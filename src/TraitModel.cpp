

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
#include "TraitBranchHistory.h"
#include "Settings.h"
#include "Log.h"
#include "Prior.h"
#include "Stat.h"

double TraitModel::mhColdness = 1.0;




TraitModel::TraitModel(MbRandom* ranptr, Tree* tp, Settings* sp, Prior* Pr)
{

    _lastLH = 0.0;

    // reduce weird autocorrelation of values at start by calling RNG a few times...
    for (int i = 0; i < 100; i++)
        ranptr->uniformRv();

    ran = ranptr;
    treePtr = tp;
    sttings = sp;
    cprior = Pr;
    
    
    treePtr->getTotalMapLength(); // total map length (required to set priors)

    // Set parameter values for model object, including priors etc.

    gen = 0;

    /*  ******************  */
    /* MCMC proposal/tuning parameters  */
    /*
     I am keeping all of these parameters with apparent double-definitions:
     They are stored both in Settings (sttings)
     but also as private data for class Model

     The rationale for this is that I may keep the class Settings as a
     one-way class that holds INITIAL settings
     But may incorporate some auto-tuning of parameters to make the MCMC
     more efficient, and thus the within-Model parameters will be updated
     but not the Settings parameters.

     This could all change, though.
     */

    _updateBetaScale = sttings->getUpdateBetaScale();
    _updateBetaShiftScale = sttings->getUpdateBetaShiftScale();

    // Node state scale is relative to the standard deviation
    // of the trait values (located in the tree terminal nodes)
    double std_dev_traits = Stat::standard_deviation(treePtr->traitValues());
    _updateNodeStateScale = sttings->getUpdateNodeStateScale() * std_dev_traits;

    // initial values

    // Event location scale is relative to the maximum root-to-tip length
    _scale = sttings->getUpdateEventLocationScale() *
        treePtr->maxRootToTipLength();

    _updateEventRateScale = sttings->getUpdateEventRateScale();
    _localGlobalMoveRatio =
        sttings->getLocalGlobalMoveRatio(); // For Poisson process
 
    _poissonRatePrior = sttings->getPoissonRatePrior();
    

    setMinMaxTraitPriors();

    eventLambda = 1 / _poissonRatePrior; // event rate, initialized to generate expected number of 1 event

    /*  *********************   */
    /* Other parameters & tracking variables*/

    acceptCount = 0;
    rejectCount = 0;
    acceptLast = -1;

    //set up event at root node:
    double startTime = 0;

  
#ifdef NEGATIVE_SHIFT_PARAM
    
    // Constrain beta shift to be zero or less than zero.
    if (sttings->getBetaShiftInit() > 0){
        std::cout << "\n\n********* ERROR ******************" << std::endl;
        std::cout << "Initial value of beta shift (betaShiftInit) cannot" << std::endl;
        std::cout << " be positive. This parameter is constrained to negative values\n\n" << std::endl;
        exit(0);
    }
#endif
    
    
    TraitBranchEvent* x =  new TraitBranchEvent((double)sttings->getBetaInit(),
            sttings->getBetaShiftInit(), treePtr->getRoot(), treePtr, ran, startTime);
    rootEvent = x;
    lastEventModified = x;

    log() << "\nRoot beta: " << x->getBetaInit() << "\t"
          << sttings->getBetaInit() << "\tShift: "
          << x->getBetaShift() << "\n";

    // set NodeEvent of root node equal to the rootEvent:
    tp->getRoot()->getTraitBranchHistory()->setNodeEvent(rootEvent);

    //initializing all branch histories to equal the root event:
    forwardSetBranchHistories(rootEvent);

    treePtr->setMeanBranchTraitRates();

    if (sttings->getLoadEventData()) {
        log() << "\nLoading model data from file: "
              << sttings->getEventDataInfile() << "\n";
        initializeModelFromEventDataFileTrait();
    }

    setCurrLnLTraits(computeLikelihoodTraits());

    log() << "\nInitial log-likelihood: " << getCurrLnLTraits() << "\n";
    if (sttings->getSampleFromPriorOnly())
        log() << "Note that you have chosen to sample from prior only.\n";

    // this parameter only set during model-jumping.
    _logQratioJump = 0.0;

}


TraitModel::~TraitModel(void)
{

    for (std::set<TraitBranchEvent*>::iterator it = eventCollection.begin();
            it != eventCollection.end(); it++)
        delete (*it);
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

    std::ifstream infile(sttings->getEventDataInfile().c_str());
    std::cout << "Initializing model from <<" << sttings->getEventDataInfile() << ">>" <<
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

            Node* x = treePtr->getNodeMRCA(species1[i].c_str(), species2[i].c_str());
            if (x  == treePtr->getRoot()) {

                // Only including this "root time setting" to replicate previous bug - Nov 1 2013
                rootEvent->setAbsoluteTime(etime[i]);
                rootEvent->setBetaInit(beta_par1[i]);
                rootEvent->setBetaShift(beta_par2[i]);

            } else {
                double deltaT = x->getTime() - etime[i];

                double newmaptime = x->getMapStart() + deltaT;

                TraitBranchEvent* newEvent = new TraitBranchEvent(beta_par1[i], beta_par2[i], x,
                        treePtr, ran, newmaptime);
                newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
                    newEvent);
                eventCollection.insert(newEvent);
                forwardSetBranchHistories(newEvent);
                treePtr->setMeanBranchTraitRates();
            }

        } else if ((species2[i] == "NA") && (species1[i] != "NA")) {

            Node* x = treePtr->getNodeByName(species1[i].c_str());

            double deltaT = x->getTime() - etime[i];
            double newmaptime = x->getMapStart() + deltaT;

            TraitBranchEvent* newEvent = new TraitBranchEvent(beta_par1[i], beta_par2[i], x,
                    treePtr, ran, newmaptime);
            newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
                newEvent);
            eventCollection.insert(newEvent);
            forwardSetBranchHistories(newEvent);
            treePtr->setMeanBranchTraitRates();

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

    double newbeta = cprior->generateBetaInitFromPrior();
    double newBetaShift = cprior->generateBetaShiftFromPrior();
    
    
#ifdef NEGATIVE_SHIFT_PARAM
    
    newBetaShift = -fabs(newBetaShift);
    double dens_term = log(2.0);
    
#else
    double dens_term = 0.0;
    
    
#endif
    
    
    _logQratioJump = 0.0;
    
    _logQratioJump += cprior->betaInitPrior(newbeta);
    _logQratioJump += dens_term + cprior->betaShiftPrior(newBetaShift);
    
    //_logQratioJump += dens_term + ran->lnExponentialPdf(sttings->getBetaInitPrior(), newbeta);
    
    // Add log(2) [see dens_term above] because this is truncated normal distribution constrained to negative values
    //_logQratioJump += dens_term + ran->lnNormalPdf((double)0.0, sttings->getBetaShiftPrior(),newBetaShift);

    // End calculations:: now create event

    TraitBranchEvent* newEvent = new TraitBranchEvent(newbeta, newBetaShift,
            treePtr->mapEventToTree(x), treePtr, ran, x);

    // add the event to the branch history.
    //  ALWAYS done after event is added to tree.
    newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
        newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    treePtr->setMeanBranchTraitRates();

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
    
    double aa = treePtr->getRoot()->getMapStart();
    double bb = treePtr->getTotalMapLength();
    double x = ran->uniformRv(aa, bb);
    
    
    /*      ********************* */
    // Sample beta and beta shift from prior:
    
    double newbeta = cprior->generateBetaInitFromPrior();
    double newBetaShift = cprior->generateBetaShiftFromPrior();
    
    
#ifdef NEGATIVE_SHIFT_PARAM
    
    newBetaShift = -fabs(newBetaShift);
    double dens_term = log(2.0);
    
#else
    double dens_term = 0.0;
    
    
#endif
    
    
    
    // End calculations:: now create event
    
    _logQratioJump = 0.0;
    
    _logQratioJump = cprior->betaInitPrior(newbeta);
    _logQratioJump += dens_term + cprior->betaShiftPrior(newBetaShift);
    
    TraitBranchEvent* newEvent = new TraitBranchEvent(newbeta, newBetaShift,
                                                      treePtr->mapEventToTree(x), treePtr, ran, x);
    
    
    
    newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
                                                                               newEvent);
    
    eventCollection.insert(newEvent);
    
    // Event is now inserted into branch history:
    //  however, branch histories must be updated.
    
    forwardSetBranchHistories(newEvent);
    treePtr->setMeanBranchTraitRates();
    
    // Addition June17 2012
    lastEventModified = newEvent;
    
    
}


// This function for adding event with beta...

void TraitModel::addEventToTreeWithSetBeta(double beta, double bshift)
{
    
    double aa = treePtr->getRoot()->getMapStart();
    double bb = treePtr->getTotalMapLength();
    double x = ran->uniformRv(aa, bb);


    // For now, the rates of speciation and extinction are set to whatever they should be based
    // on the ancestralNodeEvent
    //Node * xnode = treePtr->mapEventToTree(x);
    //double atime = treePtr->getAbsoluteTimeFromMapTime(x);
    //TraitBranchHistory * bh = xnode->getTraitBranchHistory();
    //TraitBranchEvent * be = bh->getAncestralNodeEvent();

    //double elapsed = atime - be->getAbsoluteTime();

    // End calculations:: now create event

    TraitBranchEvent* newEvent = new TraitBranchEvent(beta, bshift,
            treePtr->mapEventToTree(x), treePtr, ran, x);

    newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
        newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);
    treePtr->setMeanBranchTraitRates();

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
    for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); i++) {
        std::cout << "event " << counter++ << "\tAddress: " << (*i) << "\t";
        std::cout << (*i)->getMapTime() << "\tNode: " << (*i)->getEventNode() << std::endl <<
             std::endl;
    }


}

TraitBranchEvent* TraitModel::chooseEventAtRandom(void)
{

    int n_events = (int)eventCollection.size();
    if (n_events == 0) {
        return NULL;
        //should ultimately throw exception here.

    } else {
        int ctr = 0;
        double xx = ran->uniformRv();
        int chosen = (int)(xx * (double)n_events);

        std::set<TraitBranchEvent*>::iterator sit = eventCollection.begin();

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
        TraitBranchEvent* chosenEvent = chooseEventAtRandom();

        // corresponding node defining branch on which event occurs
        //Node* theEventNode = chosenEvent->getEventNode();

        // this is the event preceding the chosen event: histories should be set forward from here..
        TraitBranchEvent* previousEvent =
            chosenEvent->getEventNode()->getTraitBranchHistory()->getLastEvent(chosenEvent);

        // set this history variable in case move is rejected
        lastEventModified = chosenEvent;

        chosenEvent->getEventNode()->getTraitBranchHistory()->popEventOffBranchHistory(
            chosenEvent);
		
		// Get step size for move:
		double step = ran->uniformRv(0, _scale) - 0.5*_scale;
		

		
        chosenEvent->moveEventLocal(step); // move event
        chosenEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
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

    treePtr->setMeanBranchTraitRates();


}

void TraitModel::eventGlobalMove(void)
{


    if (getNumberOfEvents() > 0) {
        TraitBranchEvent* chosenEvent = chooseEventAtRandom();

        // this is the event preceding the chosen event: histories should be set forward from here..
        TraitBranchEvent* previousEvent =
            chosenEvent->getEventNode()->getTraitBranchHistory()->getLastEvent(chosenEvent);

        // private variable
        lastEventModified = chosenEvent;

        chosenEvent->getEventNode()->getTraitBranchHistory()->popEventOffBranchHistory(
            chosenEvent);
        chosenEvent->moveEventGlobal(); // move event
        chosenEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
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

    treePtr->setMeanBranchTraitRates();

}

// used to reset position of event if move is rejected

void TraitModel::revertMovedEventToPrevious(void)
{


    //double startLH = getCurrLnLTraits();

    // Get LAST EVENT from position of event to be removed:

    TraitBranchEvent* newLastEvent =
        lastEventModified->getEventNode()->getTraitBranchHistory()->getLastEvent(
            lastEventModified);

    //BranchEvent * newLastEvent = getLastEvent(lastEventModified);

    // pop event off its new position
    lastEventModified->getEventNode()->getTraitBranchHistory()->popEventOffBranchHistory(
        lastEventModified);

    // Reset nodeptr:
    // Reset mapTime:
    lastEventModified->revertOldMapPosition();

    // Now: reset forward from lastEventModified (new position)
    //  and from newLastEvent, which holds 'last' event before old position

    lastEventModified->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
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
    treePtr->setMeanBranchTraitRates();


}



// Recursively count the number of events in the branch histories
int TraitModel::countEventsInBranchHistory(Node* p)
{
    int count = 0;
    count += p->getTraitBranchHistory()->getNumberOfBranchEvents();
    if (p->getLfDesc() != NULL)
        count += countEventsInBranchHistory(p->getLfDesc());
    if (p->getRtDesc() != NULL)
        count += countEventsInBranchHistory(p->getRtDesc());

    return count;
}

/*

 Deletes an event from tree.

 */

void TraitModel::deleteEventFromTree(TraitBranchEvent* be)
{
    
    
    if (be == rootEvent) {
        std::cout << "Can't delete root event" << std::endl;
        exit(1);
    } else {
        // erase from branch history:
        Node* currNode = (be)->getEventNode();
        
        //get event downstream of i
        TraitBranchEvent* newLastEvent =
        currNode->getTraitBranchHistory()->getLastEvent(be);
        
        lastDeletedEventMapTime = (be)->getMapTime();
        
        _lastDeletedEventBetaInit = (be)->getBetaInit();
        _lastDeletedEventBetaShift = (be)->getBetaShift();
        
        /************************/
        _logQratioJump = 0.0;
        
        _logQratioJump = cprior->betaInitPrior(_lastDeletedEventBetaInit);
        _logQratioJump += cprior->betaShiftPrior(_lastDeletedEventBetaShift);
        
        currNode->getTraitBranchHistory()->popEventOffBranchHistory((be));
        
        eventCollection.erase(be);
        
        // delete from global node set
        delete (be);
        //std::cout << "deleted..." << std::endl;
        
        forwardSetBranchHistories(newLastEvent);
        
    }
    
    
    treePtr->setMeanBranchTraitRates();
    
    
}





void TraitModel::deleteRandomEventFromTree(void)
{
    
    
    //std::cout << std::endl << std::endl << "START Delete: " << std::endl;
    //printBranchHistories(treePtr->getRoot());
    
    // can only delete event if more than root node present.
    int n_events = (int)eventCollection.size();
    
    if (eventCollection.size() > 0) {
        int counter = 0;
        double xx = ran->uniformRv();
        int chosen = (int)(xx * (double)n_events);
        
        for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin();
             i != eventCollection.end(); i++) {
            if (counter++ == chosen) {
                
                // erase from branch history:
                Node* currNode = (*i)->getEventNode();
                
                //get event downstream of i
                TraitBranchEvent* newLastEvent =
                currNode->getTraitBranchHistory()->getLastEvent((*i));
                
                lastDeletedEventMapTime = (*i)->getMapTime();
                //lastDeletedEventBeta = (*i)->getBeta();
                
                _lastDeletedEventBetaInit = (*i)->getBetaInit();
                _lastDeletedEventBetaShift = (*i)->getBetaShift();
                
                
                currNode->getTraitBranchHistory()->popEventOffBranchHistory((*i));
                
                /************************/
                _logQratioJump = 0.0;
                
                _logQratioJump += cprior->betaInitPrior(_lastDeletedEventBetaInit);
                _logQratioJump += cprior->betaShiftPrior(_lastDeletedEventBetaShift);
                
                eventCollection.erase(i);
                
                // delete from global node set
                delete (*i);
                //std::cout << "deleted..." << std::endl;
                
                //std::cout << "erased ... " << std::endl;
                
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
    treePtr->setMeanBranchTraitRates();
    
}

void TraitModel::restoreLastDeletedEvent(void)
{


    // Constructor for traitEvolution model:

   
    // Use constructor for speciation and extinction

    TraitBranchEvent* newEvent = new TraitBranchEvent((double)0.0, (double)0.0,
            treePtr->mapEventToTree(lastDeletedEventMapTime), treePtr, ran,
            lastDeletedEventMapTime);

    newEvent->setBetaInit(_lastDeletedEventBetaInit);
    newEvent->setBetaShift(_lastDeletedEventBetaShift);

    // add the event to the branch history.
    //  ALWAYS done after event is added to tree.
    newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(
        newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    treePtr->setMeanBranchTraitRates();
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

    bool gain = (ran->uniformRv() <= 0.5);
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

        treePtr->setMeanBranchTraitRates();

        double likTraits = computeLikelihoodTraits();

#endif
        // Prior density:
        newLogPrior = computeLogPrior();

        logHR = log(eventLambda) - log(K + 1.0);
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
                acceptCount++;
                acceptLast = 1;
            } else {


                // Need to get rid of event that was just gained...
                //std::cout << "Invalid event config from addEventToTree - deleting." << std::endl;
                deleteEventFromTree(lastEventModified);
                treePtr->setMeanBranchTraitRates();
                rejectCount++;
                acceptLast = 0;
            }


        } else {

            // Need to get rid of event that was just gained...
            //std::cout << "Invalid event config from addEventToTree - deleting." << std::endl;
            deleteEventFromTree(lastEventModified);
            treePtr->setMeanBranchTraitRates();
            rejectCount++;
            acceptLast = 0;


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

        treePtr->setMeanBranchTraitRates();

        double likTraits = computeLikelihoodTraits();


#endif

        double qratio = 1.0; // if loss, can only be qratio of 1.0

        if (K == 1)
            qratio = 2.0;

        // Prior density:
        newLogPrior = computeLogPrior();

        logHR = log(K) - log(eventLambda);
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

            acceptCount++;
            acceptLast = 1;


        } else {

            restoreLastDeletedEvent();

            // Trait evolution rates on branches automatically updated after restoreLastDeletedEvent()
            setCurrLnLTraits(computeLikelihoodTraits());


            //std::cout << "loss reject restored, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;
            rejectCount++;
            acceptLast = 0;


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

        bool isLocalMove = (ran->uniformRv() <= localMoveProb);
        //std::cout << "is local: " << isLocalMove << std::endl;

        if (isLocalMove) {
            // Local move, with event drawn at random
            eventLocalMove();

#ifdef NO_DATA
            double likTraits = 0;
            double PropLnLik = likTraits;
#else

            treePtr->setMeanBranchTraitRates();

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

                acceptCount++;
                acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                treePtr->setMeanBranchTraitRates();

                rejectCount++;
                acceptLast = 0;

            }



        } else {
            // global move, event drawn at random
            eventGlobalMove();
            //std::cout << "successful global move" << std::endl;
#ifdef NO_DATA
            double likTraits = 0;
            double PropLnLik = likTraits;
#else

            treePtr->setMeanBranchTraitRates();

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
                acceptCount++;
                acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                treePtr->setMeanBranchTraitRates();
                rejectCount++;
                acceptLast = 0;


            }


        }

    } else {
        // consider proposal rejected (can't move nonexistent event)
        rejectCount++;
        acceptLast = 0;
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
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    TraitBranchEvent* be = rootEvent;

    if (toUpdate > 0) {
        std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
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
    double cterm = exp( _updateEventRateScale * (ran->uniformRv() - 0.5) );
    setEventRate(cterm * oldEventRate);
    
    
    double LogPriorRatio = cprior->poissonRatePrior(getEventRate());
    LogPriorRatio -= cprior->poissonRatePrior(oldEventRate);
    
    double logProposalRatio = log(cterm);
    double logHR = LogPriorRatio + logProposalRatio;
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    //std::cout << "ER " << oldEventRate << "\t" << cterm*oldEventRate << std::endl;
    
    if (acceptMove == true) {
        // continue
        acceptCount++;
        acceptLast = 1;
    } else {
        setEventRate(oldEventRate);
        rejectCount++;
        acceptLast = 0;
    }
    
    incrementGeneration();
    //std::cout << "Leaving UpdateEventRate" << std::endl;
}



void TraitModel::updateBetaMH(void)
{
    
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    
    TraitBranchEvent* be = rootEvent;
    
    
    if (toUpdate > 0) {
        std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;
        
        be = (*myIt);
    }
    
    
    double oldRate = be->getBetaInit();
    double cterm = exp( _updateBetaScale * (ran->uniformRv() - 0.5) );
    be->setBetaInit(cterm * oldRate);
    treePtr->setMeanBranchTraitRates();
    
    
    
    
    
#ifdef NO_DATA
    double PropLnLik = 0;
#else
    
    double PropLnLik = computeLikelihoodTraits();
    
#endif
    
    double LogPriorRatio = cprior->betaInitPrior(be->getBetaInit());
    LogPriorRatio -= cprior->betaInitPrior(oldRate);
    
    
    double LogProposalRatio = log(cterm);
    
    double likeRatio = PropLnLik - getCurrLnLTraits();
    
    double logHR = likeRatio + LogPriorRatio + LogProposalRatio;
    
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    //  std::cout << getGeneration() << "\tL1: " << startLH << "\tL2: " << getCurrLnLTraits() << std::endl;
    
    
    if (acceptMove == true) {
        //std::cout << "accept: " << oldRate << "\t" << be->getBetaInit() << std::endl;
        setCurrLnLTraits(PropLnLik);
        acceptCount++;
        acceptLast = 1;
        
    } else {
        
        // revert to previous state
        _lastLH = PropLnLik;
        
        
        be->setBetaInit(oldRate);
        treePtr->setMeanBranchTraitRates();
        acceptLast = 0;
        rejectCount++;
        
    }
    
    /*if (!acceptMove){
     std::cout << std::endl;
     std::cout << startLL << "\tCurr: " << getCurrLnLTraits() << "\tcalc: " << computeLikelihoodTraits() << std::endl;
     }*/
    
    incrementGeneration();
    
}


void TraitModel::updateBetaShiftMH(void)
{
    
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    
    TraitBranchEvent* be = rootEvent;
    
    
    if (toUpdate > 0) {
        std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;
        
        be = (*myIt);
    }
    
    double oldShift = be->getBetaShift();
    double newShift = oldShift + ran->normalRv((double)0.0, _updateBetaShiftScale);
    

#ifdef NEGATIVE_SHIFT_PARAM
    
    newShift = -fabs(newShift);
    
#endif
    
    be->setBetaShift(newShift);
    treePtr->setMeanBranchTraitRates();
    
#ifdef NO_DATA
    double PropLnLik = 0;
#else
    double PropLnLik = computeLikelihoodTraits();
    
#endif
    
    double LogPriorRatio = cprior->betaShiftPrior(newShift);
    LogPriorRatio -= cprior->betaShiftPrior(oldShift);
    
    
    double LogProposalRatio = 0.0;
    
    double likeRatio = PropLnLik - getCurrLnLTraits();
    
    double logHR = likeRatio + LogPriorRatio + LogProposalRatio;
    
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    if (acceptMove == true) {
        
        setCurrLnLTraits(PropLnLik);
        acceptCount++;
        acceptLast = 1;
        
    } else {
        
        // revert to previous state
        be->setBetaShift(oldShift);
        treePtr->setMeanBranchTraitRates();
        acceptLast = 0;
        rejectCount++;
        
    }
    
    
    
    incrementGeneration();
    
    
}

void TraitModel::updateNodeStateMH(void)
{

    Node* xnode = treePtr->chooseInternalNodeAtRandom();

#ifdef NO_DATA

#else
    double oldTriadLogLik = computeTriadLikelihoodTraits(xnode);

#endif

    double oldstate = xnode->getTraitValue();
    double newstate = oldstate +
        ran->uniformRv(-1.0 * _updateNodeStateScale, _updateNodeStateScale);
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
    if (newstate > sttings->getTraitPriorMax() ||
            newstate < sttings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrLnLTraits(PropLnLik);
        acceptCount++;
        acceptLast = 1;

    } else {
        xnode->setTraitValue(oldstate);

        rejectCount++;
        acceptLast = 0;
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
    double newstate = oldstate + ran->uniformRv((-1.0 *
                      sttings->getUpdateNodeStateScale()), sttings->getUpdateNodeStateScale());
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
    if (newstate > sttings->getTraitPriorMax() ||
            newstate < sttings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrLnLTraits(PropLnLik);
        acceptCount++;
        acceptLast = 1;

    } else {
        xnode->setTraitValue(oldstate);

        rejectCount++;
        acceptLast = 0;
    }

    incrementGeneration();





}


void TraitModel::updateDownstreamNodeStatesMH(Node* xnode)
{

    // Get list of internal node descendants from node * x
    // update each (or some fraction thereof).

    treePtr->setTempInternalNodeArray(xnode);
    for (int i = 0; i < 100; i++)
        updateNodeStateMH(treePtr->getRandomNodeFromTempArray());

    treePtr->clearTempNodeArray();


}


double TraitModel::computeLikelihoodTraits(void)
{

    double LnL = 0.0;

    //Node * tmpnode = treePtr->getRoot()->getLfDesc();

    if (sttings->getSampleFromPriorOnly())
        return 0.0;

#ifdef NO_DATA
    LnL = 0.0;
#else
    int numNodes = treePtr->getNumberOfNodes();

    // iterate over non-root nodes and compute LnL

    for (int i = 0; i < numNodes; i++) {
        Node* xnode = treePtr->getNodeFromDownpassSeq(i);
        if ( (xnode != treePtr->getRoot()) && (xnode->getCanHoldEvent() == true) ) {


            double var = xnode->getBrlen() * xnode->getMeanBeta();

            // change in phenotype:
            double delta = xnode->getTraitValue() - xnode->getAnc()->getTraitValue();

            LnL += ran->lnNormalPdf(0, var, delta);

            //std::cout << xnode << "dz: " << delta << "\tT: " << xnode->getBrlen() << "\tRate: " << xnode->getMeanBeta();
            //std::cout << "\tLf: " << ran->lnNormalPdf(0, var, delta) << std::endl;

            /*if (xnode == tmpnode){
                std::cout << tmpnode->getTraitBranchHistory()->getAncestralNodeEvent()->getBetaInit();
                std::cout << "\tDelta: " << delta << "\tvar: " << var << "\tLL: " << ran->lnNormalPdf(0, var, delta);
                std::cout << "\tBeta: " << xnode->getMeanBeta()  << std::endl;
            }*/
        }

    }

#endif

    return LnL;

}

double TraitModel::computeTriadLikelihoodTraits(Node* x)
{


    if (sttings->getSampleFromPriorOnly())
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
            logL += ran->lnNormalPdf(0,
                                     (x->getLfDesc()->getBrlen() * x->getLfDesc()->getMeanBeta()), delta);
        }


        if (x->getRtDesc()->getCanHoldEvent() == true) {
            // computation for right descendant branch
            double delta = x->getRtDesc()->getTraitValue() - x->getTraitValue();
            logL += ran->lnNormalPdf(0,
                                     (x->getRtDesc()->getBrlen() * x->getRtDesc()->getMeanBeta()), delta);


        }


        // computation for ancestral branch (unless == root)

        if (x != treePtr->getRoot()) {

            double delta = x->getTraitValue() - x->getAnc()->getTraitValue();
            logL += ran->lnNormalPdf(0, (x->getBrlen() * x->getMeanBeta()), delta);

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
    
    logPrior += cprior->betaInitPrior(rootEvent->getBetaInit());
    logPrior += dens_term + cprior->betaShiftPrior(rootEvent->getBetaShift());
    
    for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin();
         i != eventCollection.end(); i++) {
        
        logPrior += cprior->betaInitPrior((*i)->getBetaInit());
        logPrior += dens_term + cprior->betaShiftPrior((*i)->getBetaShift());
        
    }
    
    // and prior on number of events:
    
    logPrior += cprior->poissonRatePrior(getEventRate());
    
    return logPrior;
    
}





bool TraitModel::acceptMetropolisHastings(const double lnR)
{
    const double r = safeExponentiation(TraitModel::mhColdness * lnR);
    return (ran->uniformRv() < r);
}



void TraitModel::initializeBranchHistories(Node* x)
{
    //std::cout << x << std::endl;
    x->getTraitBranchHistory()->setNodeEvent(rootEvent);

    if (x->getAnc() != NULL)
        x->getTraitBranchHistory()->setAncestralNodeEvent(rootEvent);

    if (x->getLfDesc() != NULL)
        initializeBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        initializeBranchHistories(x->getRtDesc());


}



void TraitModel::printStartAndEndEventStatesForBranch(Node* x)
{

    if (x != treePtr->getRoot()) {
        std::cout << "Node: " << x << "\tAnc: " <<
             x->getTraitBranchHistory()->getAncestralNodeEvent();
        std::cout << "\tevent: " << x->getTraitBranchHistory()->getNodeEvent() << std::endl;
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

void TraitModel::forwardSetBranchHistories(TraitBranchEvent* x)
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




    if (x == rootEvent) {
        forwardSetHistoriesRecursive(myNode->getLfDesc());
        forwardSetHistoriesRecursive(myNode->getRtDesc());

    } else if (x == myNode->getTraitBranchHistory()->getLastEvent()) {
        // If TRUE, x is the most tip-wise event on branch.
        myNode->getTraitBranchHistory()->setNodeEvent(x);

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
    TraitBranchEvent* lastEvent =
        p->getAnc()->getTraitBranchHistory()->getNodeEvent();
    // set the ancestor equal to the event state of parent node:
    p->getTraitBranchHistory()->setAncestralNodeEvent(lastEvent);

    // if no events on the branch, go down to descendants and do same thing
    //  otherwise, process terminates (because it hits another event on branch
    if (p->getTraitBranchHistory()->getNumberOfBranchEvents() == 0) {
        p->getTraitBranchHistory()->setNodeEvent(lastEvent);
        if (p->getLfDesc() != NULL)
            forwardSetHistoriesRecursive(p->getLfDesc());
        if (p->getRtDesc() != NULL)
            forwardSetHistoriesRecursive(p->getRtDesc());
    }

}




void TraitModel::printBranchHistories(Node* x)
{

    if (x != treePtr->getRoot()) {
        std::cout << "Node: " << x;
        std::cout << "\t#Events: " << x->getTraitBranchHistory()->getNumberOfBranchEvents()
             << "\tStart: ";
        std::cout << x->getTraitBranchHistory()->getAncestralNodeEvent() << "\tEnd: ";
        std::cout << x->getTraitBranchHistory()->getNodeEvent() << std::endl;

    }
    if (x->getLfDesc() != NULL)
        printBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        printBranchHistories(x->getRtDesc());



}



double  TraitModel::getMHacceptanceRate(void)
{

    double arate = (double)acceptCount / ((double)acceptCount +
                                          (double)rejectCount);


    return arate;

}

void  TraitModel::resetMHacceptanceParameters(void)
{
    acceptCount = 0;
    rejectCount = 0;
    
}




TraitBranchEvent* TraitModel::getEventByIndex(int x)
{

    //int ctr = 0;
    std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
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
    count += (int)rootEvent->getIsEventTimeVariable();
    for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin();
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


    TraitBranchEvent* be = rootEvent;
    Node* xl = treePtr->getRoot()->getRandomLeftTipNode();
    Node* xr = treePtr->getRoot()->getRandomRightTipNode();
    ss << xl->getName() << "," << xr->getName() << "," << be->getAbsoluteTime() <<
       ",";
    ss << be->getBetaInit() << "," << be->getBetaShift();

    if (eventCollection.size() > 0) {
        for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin();
                i != eventCollection.end(); i++) {

            ss << "\n" << getGeneration() << ",";
            be = (*i);
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


bool TraitModel::isEventConfigurationValid(TraitBranchEvent* be)
{
    //std::cout << "enter isEventConfigValid" << std::endl;
    bool isValidConfig = false;

    if (be->getEventNode() == treePtr->getRoot()) {
        Node* rt = treePtr->getRoot()->getRtDesc();
        Node* lf = treePtr->getRoot()->getLfDesc();
        if (rt->getTraitBranchHistory()->getNumberOfBranchEvents() > 0 &&
                lf->getTraitBranchHistory()->getNumberOfBranchEvents() > 0) {
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

        if (anc == treePtr->getRoot())
            badsum++;
        else if (anc->getTraitBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;
        else {
            // nothing;
        }

        // test lf desc:
        if (lf->getTraitBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        // test rt desc
        if (rt->getTraitBranchHistory()->getNumberOfBranchEvents() > 0)
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

    TraitBranchEvent* be = rootEvent;
    std::cout << "RtBt: " << be->getBetaInit() << "\tSf: " << be->getBetaShift() <<
         "\tAtime:" << be->getAbsoluteTime() << std::endl;
    int ctr = 0;
    for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); i++) {
        be = (*i);
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

    int nnodes = treePtr->getNumberOfNodes();
    std::vector<double> tvec;
    for (int i = 0; i < nnodes; i++) {
        Node* xnode = treePtr->getNodeFromDownpassSeq(i);
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
    sttings->setTraitPriorMin(minprior);
    sttings->setTraitPriorMax(maxprior);

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
