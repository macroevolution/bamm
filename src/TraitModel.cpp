// defining this macro constrains analysis to NEGATIVE values
// for the beta shift parameter

//#define NEGATIVE_SHIFT_PARAM
#undef NEGATIVE_SHIFT_PARAM

#undef DEBUG  // This is a problem.


#include "TraitModel.h"


#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm>

#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "TraitBranchEvent.h"
#include "Settings.h"
#include "Log.h"
#include "Prior.h"
#include "Stat.h"


TraitModel::TraitModel(MbRandom* rng, Tree* tree, Settings* settings,
    Prior* prior) : Model(rng, tree, settings, prior)
{
    _lastLH = 0.0;

    _updateBetaScale = _settings->getUpdateBetaScale();
    _updateBetaShiftScale = _settings->getUpdateBetaShiftScale();

    // Node state scale is relative to the standard deviation
    // of the trait values (located in the tree terminal nodes)
    double sd_traits = Stat::standard_deviation(_tree->traitValues());
    _updateNodeStateScale = _settings->getUpdateNodeStateScale() * sd_traits;

    setMinMaxTraitPriors();

#ifdef NEGATIVE_SHIFT_PARAM
    // Constrain beta shift to be zero or less than zero.
    if (_settings->getBetaShiftInit() > 0) {
        log(Error) << "Initial value of beta shift (betaShiftInit) cannot be\n"
            << "positive. This parameter is constrained to negative values\n";
        std::exit(1);
    }
#endif
    
    BranchEvent* x = new TraitBranchEvent
        (_settings->getBetaInit(), _settings->getBetaShiftInit(),
            _tree->getRoot(), _tree, _rng, 0);
    _rootEvent = x;
    lastEventModified = x;

    TraitBranchEvent* traitRootEvent =
        static_cast<TraitBranchEvent*>(_rootEvent);

    log() << "\nRoot beta: " << traitRootEvent->getBetaInit() << "\t"
          << _settings->getBetaInit() << "\t"
          << "Shift: " << traitRootEvent->getBetaShift() << "\n";

    // Set NodeEvent of root node equal to the _rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event:
    forwardSetBranchHistories(_rootEvent);

    _tree->setMeanBranchTraitRates();

    if (_settings->getLoadEventData()) {
        log() << "\nLoading model data from file: "
              << _settings->getEventDataInfile() << "\n";
        initializeModelFromEventDataFileTrait();
    }

    setCurrLnLTraits(computeLikelihoodTraits());

    log() << "\nInitial log-likelihood: " << getCurrLnLTraits() << "\n";
    if (_settings->getSampleFromPriorOnly()) {
        log() << "Note that you have chosen to sample from prior only.\n";
    }

    // This parameter only set during model-jumping.
    _logQratioJump = 0.0;
}


TraitModel::~TraitModel(void)
{
    for (std::set<BranchEvent*>::iterator it = eventCollection.begin();
            it != eventCollection.end(); ++it)
        delete *it;
}

/*
 Adds event to tree based on reference map value
 -adds to branch history set
 -inserts into Model::eventCollection



 */

void TraitModel::initializeModelFromEventDataFileTrait(void)
{
    // Part 1. Read from file

    // Assumes each parameter is on a new line:
    // In order:
    // sp1, sp2, time (absolute), beta0, shiftparm
    // 5 parameters, so k * 5 lines in file, where k is number of events (k >= 1, if root included)

    std::ifstream infile(_settings->getEventDataInfile().c_str());
    std::cout << "Initializing model from <<" << _settings->getEventDataInfile() << ">>" <<
         std::endl;
    std::vector<std::string> species1;
    std::vector<std::string> species2;
    std::vector<double> etime;
    std::vector<double> beta_par1;
    std::vector<double> beta_par2;

    if (!infile.good()) {
        std::cout << "Bad Filename. Exiting\n" << std::endl;
        exit(1);
    }


    std::string tempstring;

    while (infile) {
        getline(infile, tempstring, '\t');
        species1.push_back(tempstring);

        getline(infile, tempstring, '\t');
        species2.push_back(tempstring);

        getline(infile, tempstring, '\t');
        etime.push_back(atof(tempstring.c_str()));

        getline(infile, tempstring, '\t');
        beta_par1.push_back(atof(tempstring.c_str()));

        getline(infile, tempstring);
        beta_par2.push_back(atof(tempstring.c_str()));

        if (infile.peek() == EOF)
            break;
    }

    infile.close();

    std::cout << "Read a total of " << species1.size() << " events" << std::endl;
    for (std::vector<std::string>::size_type i = 0; i < species1.size(); i++)
        std::cout << species1[i] << "\t" << species2[i] << "\t" << etime[i] << "\t" <<
             beta_par1[i] << "\t" << beta_par2[i] << std::endl;


    for (std::vector<std::string>::size_type i = 0; i < species1.size(); i++) {
        std::cout << std::endl << "MRCA of : " <<  species1[i] << "\t" << species2[i] << std::endl;
        if ((species2[i] != "NA") && (species1[i] != "NA")) {

            Node* x = _tree->getNodeMRCA(species1[i].c_str(), species2[i].c_str());
            if (x  == _tree->getRoot()) {

                // Only including this "root time setting" to replicate previous bug - Nov 1 2013
                _rootEvent->setAbsoluteTime(etime[i]);

                TraitBranchEvent* re =
                    static_cast<TraitBranchEvent*>(_rootEvent);
                re->setBetaInit(beta_par1[i]);
                re->setBetaShift(beta_par2[i]);

            } else {
                double deltaT = x->getTime() - etime[i];

                double newmaptime = x->getMapStart() + deltaT;

                TraitBranchEvent* newEvent =
                    new TraitBranchEvent(beta_par1[i], beta_par2[i], x, _tree,
                            _rng, newmaptime);
                newEvent->getEventNode()->getBranchHistory()->
                    addEventToBranchHistory( newEvent);
                eventCollection.insert(newEvent);
                forwardSetBranchHistories(newEvent);
                _tree->setMeanBranchTraitRates();
            }

        } else if ((species2[i] == "NA") && (species1[i] != "NA")) {

            Node* x = _tree->getNodeByName(species1[i].c_str());

            double deltaT = x->getTime() - etime[i];
            double newmaptime = x->getMapStart() + deltaT;

            TraitBranchEvent* newEvent = new TraitBranchEvent(beta_par1[i], beta_par2[i], x,
                    _tree, _rng, newmaptime);
            newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
                newEvent);
            eventCollection.insert(newEvent);
            forwardSetBranchHistories(newEvent);
            _tree->setMeanBranchTraitRates();

        } else {
            std::cout << "Error in Model::initializeModelFromEventDataFile" << std::endl;
            exit(1);
        }

        //std::cout << i << "\t" << computeLikelihoodBranches() << "\t" << computeLogPrior() << std::endl;
    }

    std::cout << "Added " << eventCollection.size() <<
         " pre-defined events to tree, plus root event" << std::endl;
    printEventData();
}




void TraitModel::addEventToTree(double x)
{


    // Sample beta and beta shift from prior:

    double newbeta = _prior->generateBetaInitFromPrior();
    double newBetaShift = _prior->generateBetaShiftFromPrior();
    
    
#ifdef NEGATIVE_SHIFT_PARAM
    
    newBetaShift = -fabs(newBetaShift);
    double dens_term = log(2.0);
    
#else
    double dens_term = 0.0;
    
    
#endif
    
    
    _logQratioJump = 0.0;
    
    _logQratioJump += _prior->betaInitPrior(newbeta);
    _logQratioJump += dens_term + _prior->betaShiftPrior(newBetaShift);
    
    //_logQratioJump += dens_term + _rng->lnExponentialPdf(_settings->getBetaInitPrior(), newbeta);
    
    // Add log(2) [see dens_term above] because this is truncated normal distribution constrained to negative values
    //_logQratioJump += dens_term + _rng->lnNormalPdf((double)0.0, _settings->getBetaShiftPrior(),newBetaShift);

    // End calculations:: now create event

    BranchEvent* newEvent = new TraitBranchEvent(newbeta, newBetaShift,
            _tree->mapEventToTree(x), _tree, _rng, x);

    // add the event to the branch history.
    //  ALWAYS done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
        newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    _tree->setMeanBranchTraitRates();

    // Addition June17 2012
    lastEventModified = newEvent;

}


/*
 Adds event to tree based on uniform RV
 -adds to branch history set
 -inserts into Model::eventCollection



 */


void TraitModel::addEventToTree(void)
{
    
    double aa = _tree->getRoot()->getMapStart();
    double bb = _tree->getTotalMapLength();
    double x = _rng->uniformRv(aa, bb);
    
    
    /*      ********************* */
    // Sample beta and beta shift from prior:
    
    double newbeta = _prior->generateBetaInitFromPrior();
    double newBetaShift = _prior->generateBetaShiftFromPrior();
    
    
#ifdef NEGATIVE_SHIFT_PARAM
    
    newBetaShift = -fabs(newBetaShift);
    double dens_term = log(2.0);
    
#else
    double dens_term = 0.0;
    
    
#endif
    
    
    
    // End calculations:: now create event
    
    _logQratioJump = 0.0;
    
    _logQratioJump = _prior->betaInitPrior(newbeta);
    _logQratioJump += dens_term + _prior->betaShiftPrior(newBetaShift);
    
    TraitBranchEvent* newEvent = new TraitBranchEvent(newbeta, newBetaShift,
                                                      _tree->mapEventToTree(x), _tree, _rng, x);
    
    
    
    newEvent->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(newEvent);
    
    eventCollection.insert(newEvent);
    
    // Event is now inserted into branch history:
    //  however, branch histories must be updated.
    
    forwardSetBranchHistories(newEvent);
    _tree->setMeanBranchTraitRates();
    
    // Addition June17 2012
    lastEventModified = newEvent;
    
    
}


// This function for adding event with beta...

void TraitModel::addEventToTreeWithSetBeta(double beta, double bshift)
{
    
    double aa = _tree->getRoot()->getMapStart();
    double bb = _tree->getTotalMapLength();
    double x = _rng->uniformRv(aa, bb);


    // For now, the rates of speciation and extinction are set to whatever they should be based
    // on the ancestralNodeEvent
    //Node * xnode = _tree->mapEventToTree(x);
    //double atime = _tree->getAbsoluteTimeFromMapTime(x);
    //TraitBranchHistory * bh = xnode->getTraitBranchHistory();
    //TraitBranchEvent * be = bh->getAncestralNodeEvent();

    //double elapsed = atime - be->getAbsoluteTime();

    // End calculations:: now create event

    BranchEvent* newEvent = new TraitBranchEvent(beta, bshift,
            _tree->mapEventToTree(x), _tree, _rng, x);

    newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
        newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);
    _tree->setMeanBranchTraitRates();

    // Addition June17 2012
    lastEventModified = newEvent;

}



void TraitModel::printEvents(void)
{

    // for each event:
    //  print:  maptime
    //          nodeptr
    //
    int n_events = (int)eventCollection.size();
    std::cout << "N_events: " << n_events << std::endl;
    int counter = 1;
    for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); ++i) {
        std::cout << "event " << counter++ << "\tAddress: " << (*i) << "\t";
        std::cout << (*i)->getMapTime() << "\tNode: " << (*i)->getEventNode() << std::endl <<
             std::endl;
    }


}

BranchEvent* TraitModel::chooseEventAtRandom(void)
{

    int n_events = (int)eventCollection.size();
    if (n_events == 0) {
        return NULL;
        //should ultimately throw exception here.

    } else {
        int ctr = 0;
        double xx = _rng->uniformRv();
        int chosen = (int)(xx * (double)n_events);

        std::set<BranchEvent*>::iterator sit = eventCollection.begin();

        for (int i = 0; i < chosen; i++) {
            sit++;
            ctr++;
        }
        return (*sit);
    }


}





/*
 void Model::eventLocalMove(void)

 IF events are on tree:
 choose event at random
 move locally and forward set branch histories etc.
 should also store previous event information to revert to previous

 */

void TraitModel::eventLocalMove(void)
{




    if (getNumberOfEvents() > 0) {

        // the event to be moved
        BranchEvent* chosenEvent = chooseEventAtRandom();

        // corresponding node defining branch on which event occurs
        //Node* theEventNode = chosenEvent->getEventNode();

        // this is the event preceding the chosen event: histories should be set forward from here..
        BranchEvent* previousEvent =
            chosenEvent->getEventNode()->getBranchHistory()->getLastEvent(chosenEvent);

        // set this history variable in case move is rejected
        lastEventModified = chosenEvent;

        chosenEvent->getEventNode()->getBranchHistory()->popEventOffBranchHistory(
            chosenEvent);
		
		// Get step size for move:
		double step = _rng->uniformRv(0, _scale) - 0.5*_scale;
		

		
        chosenEvent->moveEventLocal(step); // move event
        chosenEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
            chosenEvent);

        // Get last event from the theEventNode
        //      forward set its history
        //  Then go to the "moved" event and forward set its history


        //BranchEvent* newLastEvent = getLastEvent(theEventNode);


        forwardSetBranchHistories(previousEvent);
        forwardSetBranchHistories(chosenEvent);

    }
    // else no events to move
    //std::cout << "leave localMove" << std::endl;

    _tree->setMeanBranchTraitRates();


}

void TraitModel::eventGlobalMove(void)
{


    if (getNumberOfEvents() > 0) {
        BranchEvent* chosenEvent = chooseEventAtRandom();

        // this is the event preceding the chosen event: histories should be set forward from here..
        BranchEvent* previousEvent =
            chosenEvent->getEventNode()->getBranchHistory()->getLastEvent(chosenEvent);

        // private variable
        lastEventModified = chosenEvent;

        chosenEvent->getEventNode()->getBranchHistory()->popEventOffBranchHistory(
            chosenEvent);
        chosenEvent->moveEventGlobal(); // move event
        chosenEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
            chosenEvent);

        // Get last event from the theEventNode
        //      forward set its history
        //  Then go to the "moved" event and forward set its history
        //BranchEvent* newLastEvent = getLastEvent(theEventNode);
        //std::cout << "EGM2: moving " << chosenEvent << "\tLastEvent: " << previousEvent <<  "\tBLE: " << newLastEvent << std::endl;


        forwardSetBranchHistories(previousEvent);
        forwardSetBranchHistories(chosenEvent);

    }
    //std::cout << "leave globalMove" << std::endl;

    _tree->setMeanBranchTraitRates();

}

// used to reset position of event if move is rejected

void TraitModel::revertMovedEventToPrevious(void)
{


    //double startLH = getCurrLnLTraits();

    // Get LAST EVENT from position of event to be removed:

    BranchEvent* newLastEvent =
        lastEventModified->getEventNode()->getBranchHistory()->getLastEvent(
            lastEventModified);

    //BranchEvent * newLastEvent = getLastEvent(lastEventModified);

    // pop event off its new position
    lastEventModified->getEventNode()->getBranchHistory()->popEventOffBranchHistory(
        lastEventModified);

    // Reset nodeptr:
    // Reset mapTime:
    lastEventModified->revertOldMapPosition();

    // Now: reset forward from lastEventModified (new position)
    //  and from newLastEvent, which holds 'last' event before old position

    lastEventModified->getEventNode()->getBranchHistory()->addEventToBranchHistory(
        lastEventModified);

    // Forward set from new position
    forwardSetBranchHistories(newLastEvent);

    // forward set from event immediately rootwards from previous position:
    forwardSetBranchHistories(lastEventModified);



    // Set lastEventModified to NULL,
    //  because it has already been reset.
    //  Future implementations should check whether this is NULL
    //  before attempting to use it to set event

    lastEventModified = NULL;

    // Reset speciaton-extinction on branches
    _tree->setMeanBranchTraitRates();


}



// Recursively count the number of events in the branch histories
int TraitModel::countEventsInBranchHistory(Node* p)
{
    int count = 0;
    count += p->getBranchHistory()->getNumberOfBranchEvents();
    if (p->getLfDesc() != NULL)
        count += countEventsInBranchHistory(p->getLfDesc());
    if (p->getRtDesc() != NULL)
        count += countEventsInBranchHistory(p->getRtDesc());

    return count;
}

/*

 Deletes an event from tree.

 */

void TraitModel::deleteEventFromTree(BranchEvent* be)
{
    
    
    if (be == _rootEvent) {
        std::cout << "Can't delete root event" << std::endl;
        exit(1);
    } else {
        // erase from branch history:
        Node* currNode = (be)->getEventNode();
        
        //get event downstream of i
        BranchEvent* newLastEvent =
            currNode->getBranchHistory()->getLastEvent(be);
        
        _lastDeletedEventMapTime = (be)->getMapTime();
        
        TraitBranchEvent* event = static_cast<TraitBranchEvent*>(be);
        _lastDeletedEventBetaInit = event->getBetaInit();
        _lastDeletedEventBetaShift = event->getBetaShift();
        
        /************************/
        _logQratioJump = 0.0;
        
        _logQratioJump = _prior->betaInitPrior(_lastDeletedEventBetaInit);
        _logQratioJump += _prior->betaShiftPrior(_lastDeletedEventBetaShift);
        
        currNode->getBranchHistory()->popEventOffBranchHistory((be));
        
        eventCollection.erase(be);
        
        // delete from global node set
        delete (be);
        //std::cout << "deleted..." << std::endl;
        
        forwardSetBranchHistories(newLastEvent);
        
    }
    
    
    _tree->setMeanBranchTraitRates();
    
    
}





void TraitModel::deleteRandomEventFromTree(void)
{
    
    
    //std::cout << std::endl << std::endl << "START Delete: " << std::endl;
    //printBranchHistories(_tree->getRoot());
    
    // can only delete event if more than root node present.
    int n_events = (int)eventCollection.size();
    
    if (eventCollection.size() > 0) {
        int counter = 0;
        double xx = _rng->uniformRv();
        int chosen = (int)(xx * (double)n_events);
        
        for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
             i != eventCollection.end(); ++i) {
            if (counter++ == chosen) {
                
                // erase from branch history:
                Node* currNode = (*i)->getEventNode();
                
                //get event downstream of i
                BranchEvent* newLastEvent =
                    currNode->getBranchHistory()->getLastEvent((*i));
                
                _lastDeletedEventMapTime = (*i)->getMapTime();
                //lastDeletedEventBeta = (*i)->getBeta();
                
                TraitBranchEvent* be = static_cast<TraitBranchEvent*>(*i);
                _lastDeletedEventBetaInit = be->getBetaInit();
                _lastDeletedEventBetaShift = be->getBetaShift();
                
                
                currNode->getBranchHistory()->popEventOffBranchHistory((*i));
                
                /************************/
                _logQratioJump = 0.0;
                
                _logQratioJump += _prior->betaInitPrior(_lastDeletedEventBetaInit);
                _logQratioJump += _prior->betaShiftPrior(_lastDeletedEventBetaShift);
                
                delete *i;
                eventCollection.erase(i);
                
                // reset forward history from last event:
                //BranchEvent* lastEvent = getLastEvent(currNode);
                //forwardSetBranchHistories(lastEvent);
                forwardSetBranchHistories(newLastEvent);
                
                //std::cout << "forward set..." << std::endl;
                
                // is this correctly setting branch histories??
                
                break;
            }
        }
    }
    _tree->setMeanBranchTraitRates();
    
}

void TraitModel::restoreLastDeletedEvent(void)
{


    // Constructor for traitEvolution model:

   
    // Use constructor for speciation and extinction

    TraitBranchEvent* newEvent = new TraitBranchEvent((double)0.0, (double)0.0,
            _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _rng,
            _lastDeletedEventMapTime);

    newEvent->setBetaInit(_lastDeletedEventBetaInit);
    newEvent->setBetaShift(_lastDeletedEventBetaShift);

    // add the event to the branch history.
    //  ALWAYS done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
        newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    _tree->setMeanBranchTraitRates();
    //setCurrLnLTraits(computeLikelihoodTraits());

}



void TraitModel::changeNumberOfEventsMH(void)
{
    double oldLogPrior = computeLogPrior();
    double newLogPrior = 0.0;
    int currState = (int)eventCollection.size();
    int proposedState = 0;
    bool acceptMove = false;
    double oldLogLikelihood = getCurrLnLTraits();

    double logHR = 0.0;

    // current number of events on tree:
    double K = (double)currState;

    // Propose gains & losses equally if not on boundary (n = 0) events:

    bool gain = (_rng->uniformRv() <= 0.5);
    if (K == 0) {
        // set event to gain IF on boundary
        gain = true;
    }

    // now to adjust acceptance ratios:

    if (gain) {

        proposedState = currState + 1;

        double qratio = 1.0;
        if (K == 0) {
            // no events on tree
            // can only propose gains.
            qratio = 0.5;
        } else {
            // DO nothing.
        }

#ifdef NO_DATA
        double likTraits = 0.0;
        double PropLnLik = likTraits;

#else


        addEventToTree();

        _tree->setMeanBranchTraitRates();

        double likTraits = computeLikelihoodTraits();

#endif
        // Prior density:
        newLogPrior = computeLogPrior();

        logHR = log(_eventRate) - log(K + 1.0);
        logHR += log(qratio);

        double likeRatio = (likTraits - oldLogLikelihood);
        logHR += likeRatio;

        logHR += (newLogPrior - oldLogPrior);

        logHR -= _logQratioJump;


        if (std::isinf(likeRatio) ) {

        } else
            acceptMove = acceptMetropolisHastings(logHR);



        if (acceptMove) {

            bool isValidConfig = isEventConfigurationValid(lastEventModified);

            if (isValidConfig) {

                setCurrLnLTraits(likTraits);
                // Update accept/reject statistics
                _acceptCount++;
                _acceptLast = 1;
            } else {


                // Need to get rid of event that was just gained...
                //std::cout << "Invalid event config from addEventToTree - deleting." << std::endl;
                deleteEventFromTree(lastEventModified);
                _tree->setMeanBranchTraitRates();
                _rejectCount++;
                _acceptLast = 0;
            }


        } else {

            // Need to get rid of event that was just gained...
            //std::cout << "Invalid event config from addEventToTree - deleting." << std::endl;
            deleteEventFromTree(lastEventModified);
            _tree->setMeanBranchTraitRates();
            _rejectCount++;
            _acceptLast = 0;


        }
        //setCurrLnLTraits(computeLikelihoodTraits());




    } else {

        // LOSS of event:

        deleteRandomEventFromTree();

        proposedState = currState - 1;

        //std::cout << "loss: LH after deleteRandomEvent" << computeLikelihoodBranches() << std::endl;

#ifdef NO_DATA
        double likTraits = 0.0;
        double PropLnLik = likTraits;

#else

        _tree->setMeanBranchTraitRates();

        double likTraits = computeLikelihoodTraits();


#endif

        double qratio = 1.0; // if loss, can only be qratio of 1.0

        if (K == 1)
            qratio = 2.0;

        // Prior density:
        newLogPrior = computeLogPrior();

        logHR = log(K) - log(_eventRate);
        logHR += log(qratio);

        double likeRatio = (likTraits - oldLogLikelihood);
        logHR += likeRatio;

        logHR += (newLogPrior - oldLogPrior);

        logHR += _logQratioJump;

        if (std::isinf(likeRatio) ) {

        } else
            acceptMove = acceptMetropolisHastings(logHR);


        if (acceptMove) {

            //std::cout << "loss accept, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;
            setCurrLnLTraits(likTraits);

            _acceptCount++;
            _acceptLast = 1;


        } else {

            restoreLastDeletedEvent();

            // Trait evolution rates on branches automatically updated after restoreLastDeletedEvent()
            setCurrLnLTraits(computeLikelihoodTraits());


            //std::cout << "loss reject restored, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;
            _rejectCount++;
            _acceptLast = 0;


        }
        //setCurrLnLTraits(computeLikelihoodTraits());
    }

#ifdef DEBUG

    std::cout << "gain: " << gain << "\tOP: " << oldLogPrior << "\tNP: " << newLogPrior
         << "\tLQ: " << _logQratioJump;
    std::cout << "\tAc: " << acceptMove << "\tlogHR " << logHR;
    std::cout << "ps/k:  " << K << "\t" << proposedState << std::endl;

#endif

    incrementGeneration();

}

void TraitModel::moveEventMH(void)
{


    if (eventCollection.size() > 0) {

        double localMoveProb = _localGlobalMoveRatio / (1 + _localGlobalMoveRatio);

        bool isLocalMove = (_rng->uniformRv() <= localMoveProb);
        //std::cout << "is local: " << isLocalMove << std::endl;

        if (isLocalMove) {
            // Local move, with event drawn at random
            eventLocalMove();

#ifdef NO_DATA
            double likTraits = 0;
            double PropLnLik = likTraits;
#else

            _tree->setMeanBranchTraitRates();

            double likTraits = computeLikelihoodTraits();
            double PropLnLik = likTraits;



#endif

            double likeRatio = PropLnLik - getCurrLnLTraits();
            double logHR = likeRatio;

            // No longer protecting this as const bool

            bool acceptMove = false;
            bool isValid = false;
            //std::cout << "calling isValid from moveEventMH::local" << std::endl;
            isValid = isEventConfigurationValid(lastEventModified);

            if (std::isinf(likeRatio) ) {

            } else if (isValid)
                acceptMove = acceptMetropolisHastings(logHR);
            else {
                //std::cout << "Invalid event configuration from LocalMove" << std::endl;
            }

            //const bool acceptMove = acceptMetropolisHastings(logHR);

            if (acceptMove == true) {
                setCurrLnLTraits(likTraits);

                _acceptCount++;
                _acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                _tree->setMeanBranchTraitRates();

                _rejectCount++;
                _acceptLast = 0;

            }



        } else {
            // global move, event drawn at random
            eventGlobalMove();
            //std::cout << "successful global move" << std::endl;
#ifdef NO_DATA
            double likTraits = 0;
            double PropLnLik = likTraits;
#else

            _tree->setMeanBranchTraitRates();

            double likTraits = computeLikelihoodTraits();
            double PropLnLik = likTraits;

#endif

            double likeRatio = PropLnLik - getCurrLnLTraits();
            double logHR = likeRatio;

            //const bool acceptMove = acceptMetropolisHastings(logHR);

            bool acceptMove = false;
            bool isValid = false;
            //std::cout << "calling isValid from moveEventMH::global" << std::endl;
            isValid = isEventConfigurationValid(lastEventModified);

            if (std::isinf(likeRatio) ) {

            } else if (isValid)
                acceptMove = acceptMetropolisHastings(logHR);
            else {
                //std::cout << "Invalid event configuration from GlobalMove" << std::endl;
            }

            if (acceptMove == true) {
                setCurrLnLTraits(likTraits);
                _acceptCount++;
                _acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                _tree->setMeanBranchTraitRates();
                _rejectCount++;
                _acceptLast = 0;


            }


        }

    } else {
        // consider proposal rejected (can't move nonexistent event)
        _rejectCount++;
        _acceptLast = 0;
    }


    incrementGeneration();

}




/* June 12 2012
 Select an event at random.
 If partition is time-constant
 flip state to time-variable
 If partition is time-variable
 flip state to time-constant

 */
void TraitModel::updateTimeVariablePartitionsMH(void)
{

    //int n_events = eventCollection.size() + 1;
    int toUpdate = _rng->sampleInteger(0, (int)eventCollection.size());
    BranchEvent* be = _rootEvent;

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = (*myIt);
    } else {
        // event remains as root event-
    }

    if (be->getIsEventTimeVariable()) {




    } else if (!be->getIsEventTimeVariable()) {



    } else {
        // Should not be able to get here:
        std::cout << "Invalid _isEventTimeVariable in Model::UpdateTimeVariablePartitionsMH"
             << std::endl;
        throw;
    }


}



/*

 Metropolis-Hastings step to update Poisson event rate.
 Note that changing this rate does not affect the likelihood,
 so the priors and qratio determine acceptance rate.

 */
void TraitModel::updateEventRateMH(void)
{
    //std::cout << "Entering update event rate" << std::endl;
    
    double oldEventRate = getEventRate();
    double cterm = exp( _updateEventRateScale * (_rng->uniformRv() - 0.5) );
    setEventRate(cterm * oldEventRate);
    
    
    double LogPriorRatio = _prior->poissonRatePrior(getEventRate());
    LogPriorRatio -= _prior->poissonRatePrior(oldEventRate);
    
    double logProposalRatio = log(cterm);
    double logHR = LogPriorRatio + logProposalRatio;
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    //std::cout << "ER " << oldEventRate << "\t" << cterm*oldEventRate << std::endl;
    
    if (acceptMove == true) {
        // continue
        _acceptCount++;
        _acceptLast = 1;
    } else {
        setEventRate(oldEventRate);
        _rejectCount++;
        _acceptLast = 0;
    }
    
    incrementGeneration();
    //std::cout << "Leaving UpdateEventRate" << std::endl;
}



void TraitModel::updateBetaMH(void)
{
    
    int toUpdate = _rng->sampleInteger(0, (int)eventCollection.size());
    
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);
    
    
    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;
        
        be = static_cast<TraitBranchEvent*>(*myIt);
    }
    
    
    double oldRate = be->getBetaInit();
    double cterm = exp( _updateBetaScale * (_rng->uniformRv() - 0.5) );
    be->setBetaInit(cterm * oldRate);
    _tree->setMeanBranchTraitRates();
    
    
    
    
    
#ifdef NO_DATA
    double PropLnLik = 0;
#else
    
    double PropLnLik = computeLikelihoodTraits();
    
#endif
    
    double LogPriorRatio = _prior->betaInitPrior(be->getBetaInit());
    LogPriorRatio -= _prior->betaInitPrior(oldRate);
    
    
    double LogProposalRatio = log(cterm);
    
    double likeRatio = PropLnLik - getCurrLnLTraits();
    
    double logHR = likeRatio + LogPriorRatio + LogProposalRatio;
    
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    //  std::cout << getGeneration() << "\tL1: " << startLH << "\tL2: " << getCurrLnLTraits() << std::endl;
    
    
    if (acceptMove == true) {
        //std::cout << "accept: " << oldRate << "\t" << be->getBetaInit() << std::endl;
        setCurrLnLTraits(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
        
    } else {
        
        // revert to previous state
        _lastLH = PropLnLik;
        
        
        be->setBetaInit(oldRate);
        _tree->setMeanBranchTraitRates();
        _acceptLast = 0;
        _rejectCount++;
        
    }
    
    /*if (!acceptMove){
     std::cout << std::endl;
     std::cout << startLL << "\tCurr: " << getCurrLnLTraits() << "\tcalc: " << computeLikelihoodTraits() << std::endl;
     }*/
    
    incrementGeneration();
    
}


void TraitModel::updateBetaShiftMH(void)
{
    
    int toUpdate = _rng->sampleInteger(0, (int)eventCollection.size());
    
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);
    
    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;
        
        be = static_cast<TraitBranchEvent*>(*myIt);
    }
    
    double oldShift = be->getBetaShift();
    double newShift = oldShift + _rng->normalRv((double)0.0, _updateBetaShiftScale);
    

#ifdef NEGATIVE_SHIFT_PARAM
    
    newShift = -fabs(newShift);
    
#endif
    
    be->setBetaShift(newShift);
    _tree->setMeanBranchTraitRates();
    
#ifdef NO_DATA
    double PropLnLik = 0;
#else
    double PropLnLik = computeLikelihoodTraits();
    
#endif
    
    double LogPriorRatio = _prior->betaShiftPrior(newShift);
    LogPriorRatio -= _prior->betaShiftPrior(oldShift);
    
    
    double LogProposalRatio = 0.0;
    
    double likeRatio = PropLnLik - getCurrLnLTraits();
    
    double logHR = likeRatio + LogPriorRatio + LogProposalRatio;
    
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    if (acceptMove == true) {
        
        setCurrLnLTraits(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
        
    } else {
        
        // revert to previous state
        be->setBetaShift(oldShift);
        _tree->setMeanBranchTraitRates();
        _acceptLast = 0;
        _rejectCount++;
        
    }
    
    
    
    incrementGeneration();
    
    
}

void TraitModel::updateNodeStateMH(void)
{

    Node* xnode = _tree->chooseInternalNodeAtRandom();

#ifdef NO_DATA

#else
    double oldTriadLogLik = computeTriadLikelihoodTraits(xnode);

#endif

    double oldstate = xnode->getTraitValue();
    double newstate = oldstate +
        _rng->uniformRv(-1.0 * _updateNodeStateScale, _updateNodeStateScale);
    xnode->setTraitValue(newstate);

#ifdef NO_DATA
    double PropLnLik = 0.0;
#else
    // Compute triad likelihood
    double newTriadLoglik = computeTriadLikelihoodTraits(xnode);
    double PropLnLik = getCurrLnLTraits() - oldTriadLogLik + newTriadLoglik;

#endif

    // set to zero for now...flat prior, unifor (see below).
    double LogPriorRatio = 0.0;


    double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0

    double likeRatio = PropLnLik - getCurrLnLTraits();


    double logHR = LogPriorRatio + logProposalRatio + likeRatio;
    bool acceptMove = acceptMetropolisHastings(logHR);

    // Here we do prior calculation to avoid computing infinity...
    if (newstate > _settings->getTraitPriorMax() ||
            newstate < _settings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrLnLTraits(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {
        xnode->setTraitValue(oldstate);

        _rejectCount++;
        _acceptLast = 0;
    }

    incrementGeneration();


}


void TraitModel::updateNodeStateMH(Node* xnode)
{

#ifdef NO_DATA

#else
    double oldTriadLogLik = computeTriadLikelihoodTraits(xnode);

#endif

    double oldstate = xnode->getTraitValue();
    double newstate = oldstate + _rng->uniformRv((-1.0 *
                      _settings->getUpdateNodeStateScale()), _settings->getUpdateNodeStateScale());
    xnode->setTraitValue(newstate);

#ifdef NO_DATA
    double PropLnLik = 0.0;
#else
    // Compute triad likelihood
    double newTriadLoglik = computeTriadLikelihoodTraits(xnode);
    double PropLnLik = getCurrLnLTraits() - oldTriadLogLik + newTriadLoglik;



#endif

    // set to zero for now...flat prior, unifor (see below).
    double LogPriorRatio = 0.0;


    double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0

    double likeRatio = PropLnLik - getCurrLnLTraits();


    double logHR = LogPriorRatio + logProposalRatio + likeRatio;
    bool acceptMove = acceptMetropolisHastings(logHR);

    // Here we do prior calculation to avoid computing infinity...
    if (newstate > _settings->getTraitPriorMax() ||
            newstate < _settings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrLnLTraits(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {
        xnode->setTraitValue(oldstate);

        _rejectCount++;
        _acceptLast = 0;
    }

    incrementGeneration();





}


void TraitModel::updateDownstreamNodeStatesMH(Node* xnode)
{

    // Get list of internal node descendants from node * x
    // update each (or some fraction thereof).

    _tree->setTempInternalNodeArray(xnode);
    for (int i = 0; i < 100; i++)
        updateNodeStateMH(_tree->getRandomNodeFromTempArray());

    _tree->clearTempNodeArray();


}


double TraitModel::computeLikelihoodTraits(void)
{

    double LnL = 0.0;

    //Node * tmpnode = _tree->getRoot()->getLfDesc();

    if (_settings->getSampleFromPriorOnly())
        return 0.0;

#ifdef NO_DATA
    LnL = 0.0;
#else
    int numNodes = _tree->getNumberOfNodes();

    // iterate over non-root nodes and compute LnL

    for (int i = 0; i < numNodes; i++) {
        Node* xnode = _tree->getNodeFromDownpassSeq(i);
        if ( (xnode != _tree->getRoot()) && (xnode->getCanHoldEvent() == true) ) {


            double var = xnode->getBrlen() * xnode->getMeanBeta();

            // change in phenotype:
            double delta = xnode->getTraitValue() - xnode->getAnc()->getTraitValue();

            LnL += _rng->lnNormalPdf(0, var, delta);

            //std::cout << xnode << "dz: " << delta << "\tT: " << xnode->getBrlen() << "\tRate: " << xnode->getMeanBeta();
            //std::cout << "\tLf: " << _rng->lnNormalPdf(0, var, delta) << std::endl;

            /*if (xnode == tmpnode){
                std::cout << tmpnode->getTraitBranchHistory()->getAncestralNodeEvent()->getBetaInit();
                std::cout << "\tDelta: " << delta << "\tvar: " << var << "\tLL: " << _rng->lnNormalPdf(0, var, delta);
                std::cout << "\tBeta: " << xnode->getMeanBeta()  << std::endl;
            }*/
        }

    }

#endif

    return LnL;

}

double TraitModel::computeTriadLikelihoodTraits(Node* x)
{


    if (_settings->getSampleFromPriorOnly())
        return 0.0;

#ifdef DEBUG
    std::cout << "Enter computeTriadLikelihood: Node : " << x << std::endl;
#endif

    double logL = 0.0;

    // Can only use this likelihood if node contributes to
    // likelihood of observed data

    if (x->getCanHoldEvent() == true) {

        // computation for left descendant branch:

        if (x->getLfDesc()->getCanHoldEvent() == true) {
            double delta = x->getLfDesc()->getTraitValue() - x->getTraitValue();
            logL += _rng->lnNormalPdf(0,
                                     (x->getLfDesc()->getBrlen() * x->getLfDesc()->getMeanBeta()), delta);
        }


        if (x->getRtDesc()->getCanHoldEvent() == true) {
            // computation for right descendant branch
            double delta = x->getRtDesc()->getTraitValue() - x->getTraitValue();
            logL += _rng->lnNormalPdf(0,
                                     (x->getRtDesc()->getBrlen() * x->getRtDesc()->getMeanBeta()), delta);


        }


        // computation for ancestral branch (unless == root)

        if (x != _tree->getRoot()) {

            double delta = x->getTraitValue() - x->getAnc()->getTraitValue();
            logL += _rng->lnNormalPdf(0, (x->getBrlen() * x->getMeanBeta()), delta);

        }



    }

#ifdef DEBUG
    std::cout << "Leaving computeTriadLikelihood: Node : " << x << std::endl;
#endif

    return logL;

}





double TraitModel::computeLogPrior(void)
{
    
    
    
    
#ifdef NEGATIVE_SHIFT_PARAM
    
    double dens_term = log(2.0);
    
#else
    double dens_term = 0.0;
    
#endif
    
    
    double logPrior = 0.0;

    TraitBranchEvent* re = static_cast<TraitBranchEvent*>(_rootEvent);
    
    logPrior += _prior->betaInitPrior(re->getBetaInit());
    logPrior += dens_term + _prior->betaShiftPrior(re->getBetaShift());
    
    for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
         i != eventCollection.end(); ++i) {

        TraitBranchEvent* event = static_cast<TraitBranchEvent*>(*i);
        
        logPrior += _prior->betaInitPrior(event->getBetaInit());
        logPrior += dens_term + _prior->betaShiftPrior(event->getBetaShift());
        
    }
    
    // and prior on number of events:
    
    logPrior += _prior->poissonRatePrior(getEventRate());
    
    return logPrior;
    
}





bool TraitModel::acceptMetropolisHastings(const double lnR)
{
    const double r = safeExponentiation(TraitModel::mhColdness * lnR);
    return (_rng->uniformRv() < r);
}



void TraitModel::initializeBranchHistories(Node* x)
{
    //std::cout << x << std::endl;
    x->getBranchHistory()->setNodeEvent(_rootEvent);

    if (x->getAnc() != NULL)
        x->getBranchHistory()->setAncestralNodeEvent(_rootEvent);

    if (x->getLfDesc() != NULL)
        initializeBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        initializeBranchHistories(x->getRtDesc());


}



void TraitModel::printStartAndEndEventStatesForBranch(Node* x)
{

    if (x != _tree->getRoot()) {
        std::cout << "Node: " << x << "\tAnc: " <<
             x->getBranchHistory()->getAncestralNodeEvent();
        std::cout << "\tevent: " << x->getBranchHistory()->getNodeEvent() << std::endl;
    }

    if (x->getLfDesc() != NULL)
        printStartAndEndEventStatesForBranch(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        printStartAndEndEventStatesForBranch(x->getRtDesc());
}




/*
 If this works correctly, this will take care of the following:
 1. if a new event is created or added to tree,
 this will forward set all branch histories from the insertion point

 2. If an event is deleted, you find the next event rootwards,
 and call forwardSetBranchHistories from that point. It will replace
 settings due to the deleted node with the next rootwards node.

 */

void TraitModel::forwardSetBranchHistories(BranchEvent* x)
{
    // If there is another event occurring more recent (closer to tips)
    //  do nothing. Even just sits in BranchHistory but doesn't affect
    //  state of any other nodes.


    // this seems circular, but what else to do?
    //  given an event (which references the node defining the branch on which event occurs)
    //   you get the corresponding branch history and the last event
    //   since the events will have been inserted in the correct order.

    Node* myNode = x->getEventNode();
    //std::cout << "Node: " << myNode << std::endl;

    //std::cout << std::endl << std::endl;
    //std::cout << "event in forwardSet: " << x << std::endl;

    //printEventData();




    if (x == _rootEvent) {
        forwardSetHistoriesRecursive(myNode->getLfDesc());
        forwardSetHistoriesRecursive(myNode->getRtDesc());

    } else if (x == myNode->getBranchHistory()->getLastEvent()) {
        // If TRUE, x is the most tip-wise event on branch.
        myNode->getBranchHistory()->setNodeEvent(x);

        // if myNode is not a tip:
        if (myNode->getLfDesc() != NULL && myNode->getRtDesc() != NULL) {
            forwardSetHistoriesRecursive(myNode->getLfDesc());
            forwardSetHistoriesRecursive(myNode->getRtDesc());
        }
        // else: node is a tip : do nothing.


    }
    //else: there is another more tipwise event on same branch; do nothing


}


void TraitModel::forwardSetHistoriesRecursive(Node* p)
{

    // Get event that characterizes parent node
    BranchEvent* lastEvent =
        p->getAnc()->getBranchHistory()->getNodeEvent();
    // set the ancestor equal to the event state of parent node:
    p->getBranchHistory()->setAncestralNodeEvent(lastEvent);

    // if no events on the branch, go down to descendants and do same thing
    //  otherwise, process terminates (because it hits another event on branch
    if (p->getBranchHistory()->getNumberOfBranchEvents() == 0) {
        p->getBranchHistory()->setNodeEvent(lastEvent);
        if (p->getLfDesc() != NULL)
            forwardSetHistoriesRecursive(p->getLfDesc());
        if (p->getRtDesc() != NULL)
            forwardSetHistoriesRecursive(p->getRtDesc());
    }

}




void TraitModel::printBranchHistories(Node* x)
{

    if (x != _tree->getRoot()) {
        std::cout << "Node: " << x;
        std::cout << "\t#Events: " << x->getBranchHistory()->getNumberOfBranchEvents()
             << "\tStart: ";
        std::cout << x->getBranchHistory()->getAncestralNodeEvent() << "\tEnd: ";
        std::cout << x->getBranchHistory()->getNodeEvent() << std::endl;

    }
    if (x->getLfDesc() != NULL)
        printBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        printBranchHistories(x->getRtDesc());



}



double  TraitModel::getMHacceptanceRate(void)
{

    double arate = (double)_acceptCount / ((double)_acceptCount +
                                          (double)_rejectCount);


    return arate;

}

void  TraitModel::resetMHacceptanceParameters(void)
{
    _acceptCount = 0;
    _rejectCount = 0;
    
}




BranchEvent* TraitModel::getEventByIndex(int x)
{

    //int ctr = 0;
    std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
    for (int i = 0; i <= x; i++)
        myIt++;

    return (*myIt);
}



/*
 Model::countTimeVaryingRatePartitions

 -counts number of time-varying rate partitions

 */
int TraitModel::countTimeVaryingRatePartitions(void)
{

    int count = 0;
    count += (int)_rootEvent->getIsEventTimeVariable();
    for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); i++)
        count += (int)(*i)->getIsEventTimeVariable();
    return count;
}


/*
 Write event data to file for all events "on" tree
 at a given point in the MCMC chain


 */

void TraitModel::getEventDataString(std::stringstream& ss)
{


    ss << getGeneration() << ",";


    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);
    Node* xl = _tree->getRoot()->getRandomLeftTipNode();
    Node* xr = _tree->getRoot()->getRandomRightTipNode();
    ss << xl->getName() << "," << xr->getName() << "," << be->getAbsoluteTime() <<
       ",";
    ss << be->getBetaInit() << "," << be->getBetaShift();

    if (eventCollection.size() > 0) {
        for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
                i != eventCollection.end(); ++i) {

            ss << "\n" << getGeneration() << ",";
            be = static_cast<TraitBranchEvent*>(*i);
            if (be->getEventNode()->getLfDesc() == NULL)
                ss << be->getEventNode()->getName() << ",NA,";
            else {
                Node* xl = be->getEventNode()->getRandomLeftTipNode();
                Node* xr = be->getEventNode()->getRandomRightTipNode();

                ss << xl->getName() << "," << xr->getName() << ",";
            }
            ss << be->getAbsoluteTime() << "," << be->getBetaInit() << "," <<
               be->getBetaShift();

        }



    }

}


bool TraitModel::isEventConfigurationValid(BranchEvent* be)
{
    //std::cout << "enter isEventConfigValid" << std::endl;
    bool isValidConfig = false;

    if (be->getEventNode() == _tree->getRoot()) {
        Node* rt = _tree->getRoot()->getRtDesc();
        Node* lf = _tree->getRoot()->getLfDesc();
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0 &&
                lf->getBranchHistory()->getNumberOfBranchEvents() > 0) {
            // events on both descendants of root. This fails.
            isValidConfig = false;
        } else
            isValidConfig = true;

    } else {
        int badsum = 0;

        Node* anc = be->getEventNode()->getAnc();
        Node* lf = anc->getLfDesc();
        Node* rt = anc->getRtDesc();

        //std::cout << "a: " << anc << "\tb: " << lf << "\tc: " << rt << std::endl;

        // test ancestor for events on branch:

        if (anc == _tree->getRoot())
            badsum++;
        else if (anc->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;
        else {
            // nothing;
        }

        // test lf desc:
        if (lf->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        // test rt desc
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        if (badsum == 3)
            isValidConfig = false;
        else if (badsum < 3)
            isValidConfig = true;
        else {
            std::cout << "problem in Model::isEventConfigurationValid" << std::endl;
            exit(1);
        }


    }


    //std::cout << "leaving isEventConfigValid. Value: " << isValidConfig << std::endl;
    return isValidConfig;
}



void TraitModel::printEventData(void)
{

    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);
    std::cout << "RtBt: " << be->getBetaInit() << "\tSf: " << be->getBetaShift() <<
         "\tAtime:" << be->getAbsoluteTime() << std::endl;
    int ctr = 0;
    for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); i++) {
        be = static_cast<TraitBranchEvent*>(*i);
        std::cout << ctr++ << "\tBt: " << be->getBetaInit() << "\tSt: " << be->getBetaShift()
             << "\tMap: " << be->getMapTime();
        std::cout << "\tAtime:" << be->getAbsoluteTime() << std::endl;

    }
    std::cout << std::endl;
}

/*
void TraitModel::initializeTraitParamsForNodes(void){

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++){


    }


}



*/


void TraitModel::setMinMaxTraitPriors(void)
{

    int nnodes = _tree->getNumberOfNodes();
    std::vector<double> tvec;
    for (int i = 0; i < nnodes; i++) {
        Node* xnode = _tree->getNodeFromDownpassSeq(i);
        if (xnode->getTraitValue() != 0)
            tvec.push_back(xnode->getTraitValue());
    }

    std::sort(tvec.begin(), tvec.end());

    // std::cout << "Min: " << tvec[0] << "\tMax: " << tvec[(tvec.size() - 1)] << std::endl;

    // Default here will be to use observed range +/- 20%
    double rg = tvec[(tvec.size() - 1)] - tvec[0];
    double minprior = tvec[0] - (0.2 * rg);
    double maxprior = tvec[(tvec.size() - 1)] + (0.2 * rg);

    log() << "\nMin and max phenotype limits set using observed data: \n"
          << "\t\tMin: " << minprior << "\tMax: " << maxprior << "\n";
    _settings->setTraitPriorMin(minprior);
    _settings->setTraitPriorMax(maxprior);

}


double TraitModel::safeExponentiation(double x)
{
    if (x > 0.0)
        return 1.0;
    else if (x < -100.0)
        return 0.0;
    else
        return exp(x);
}
