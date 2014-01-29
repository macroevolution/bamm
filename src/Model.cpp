/*
 *  Model.cpp
 *  rateshift
 *
 *  Created by Dan Rabosky on 12/1/11.
 *
 */

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>

#include "Model.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "BranchHistory.h"
#include "Settings.h"
#include "Log.h"

#include "Prior.h"

#define NO_DATA
#undef NO_DATA

#define RECURSIVE_NODES
#undef RECURSIVE_NODES

#define DEBUG
#undef DEBUG

#define DEBUG_TIME_VARIABLE

#define ADAPTIVE_MCMC_PROPOSAL
#undef ADAPTIVE_MCMC_PROPOSAL

#define NO_TRAIT
//#undef NO_TRAIT

// Maximum value of extinction probability on branch that will be tolerated:
//  to avoid numerical overflow issues (especially rounding to 1)
#define MAX_E_PROB 0.999
#define JUMP_VARIANCE_NORMAL 0.05


double Model::mhColdness = 1.0;

Model::Model(MbRandom* ranptr, Tree* tp, Settings* sp, Prior* pr)
{
    // reduce weird autocorrelation of values at start by calling RNG a few times...
    for (int i = 0; i < 100; i++)
        ranptr->uniformRv();

    ran = ranptr;
    treePtr = tp;
    sttings = sp;

    cprior = pr;
    
    
    
    // Initialize prior object here
    
    


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

    _updateLambdaInitScale = sttings->getUpdateLambdaInitScale();
    _updateMuInitScale = sttings->getUpdateMuInitScale();

    _updateLambdaShiftScale =
        sttings->getUpdateLambdaShiftScale(); // for normal proposal
    _updateMuShiftScale = sttings->getUpdateMuShiftScale();

    // initial values
    _lambdaInit0 = sttings->getLambdaInit0();
    _lambdaShift0 = sttings->getLambdaShift0();
    _muInit0 = sttings->getMuInit0();
    _muShift0 = sttings->getMuShift0();

    // Event location scale is relative to the maximum root-to-tip length
    _scale = sttings->getUpdateEventLocationScale() *
        treePtr->maxRootToTipLength();

    std::cout << "max root to tip length: " << treePtr->maxRootToTipLength()
              << "   _scale: " << _scale << "\n";

    _updateEventRateScale = sttings->getUpdateEventRateScale();
    _localGlobalMoveRatio =
        sttings->getLocalGlobalMoveRatio();
    
    // For Poisson process
    
    _poissonRatePrior = sttings->getPoissonRatePrior();

    _lambdaInitPrior = sttings->getLambdaInitPrior();
    _lambdaShiftPrior = sttings->getLambdaShiftPrior();


    _muInitPrior = sttings->getMuInitPrior();
    _muShiftPrior = sttings->getMuShiftPrior();

    eventLambda = 1 / sttings->getPoissonRatePrior(); // event rate, initialized to generate expected number of _targetNumber events

    //Parameter for splitting branch into pieces for numerical computation
    
    _segLength = sttings->getSegLength() * treePtr->maxRootToTipLength();
   
    std::cout << "relative segLength: " << sttings->getSegLength();
    std::cout << "   _segLength: " << _segLength << "\n";
    /*  *********************   */
    /* Other parameters & tracking variables*/

    acceptCount = 0;
    rejectCount = 0;
    acceptLast = -1;

    //set up event at root node:
    double startTime = 0;

    BranchEvent* x =  new BranchEvent(_lambdaInit0, sttings->getLambdaShift0(),
                                      _muInit0, sttings->getMuShift0(), treePtr->getRoot(), treePtr, ran, startTime);
    rootEvent = x;

    lastEventModified = x;
    // set NodeEvent of root node equal to the rootEvent:
    tp->getRoot()->getBranchHistory()->setNodeEvent(rootEvent);

    //std::cout << "Root time, after initializing B: " << getRootEvent()->getAbsoluteTime() << std::endl;

    //initializing all branch histories to equal the root event:
    forwardSetBranchHistories(rootEvent);

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

    //std::cout << "Root time, after initializing C: " << getRootEvent()->getAbsoluteTime() << std::endl;

    // 9.14.2012: initialize by previous event histories.
    if (sttings->getLoadEventData()) {
        log() << "\nLoading model data from file.\n";
        initializeModelFromEventDataFile();
    }
	
    setCurrLnLBranches(computeLikelihoodBranches());

    setCurrLnLTraits(0.0);
    log() << "\nInitial log-likelihood: " << getCurrLnLBranches() << "\n";
    if (sttings->getSampleFromPriorOnly())
        log() << "Note that you have chosen to sample from prior only.\n";

    // this parameter only set during model-jumping.
    _logQratioJump = 0.0;


}


Model::~Model(void)
{

    for (std::set<BranchEvent*>::iterator it = eventCollection.begin();
            it != eventCollection.end(); it++)
        delete (*it);
}


void Model::initializeModelFromEventDataFile(void)
{


    // Part 1. Read from file

    std::ifstream infile(sttings->getEventDataInfile().c_str());
    std::cout << "Initializing model from <<" << sttings->getEventDataInfile() << ">>" <<
         std::endl;
    std::vector<std::string> species1;
    std::vector<std::string> species2;
    std::vector<double> etime;
    std::vector<double> lam_par1;
    std::vector<double> lam_par2;
    std::vector<double> mu_par1;
    std::vector<double> mu_par2;

    if (!infile.good()) {
        std::cout << "Bad Filename. Exiting\n" << std::endl;
        exit(1);
    }


    std::string tempstring;

    while (infile) {
        //string tempstring;
        getline(infile, tempstring, '\t');
        //std::cout << tempstring.c_str() << "\t";
        species1.push_back(tempstring);

        getline(infile, tempstring, '\t');
        //std::cout << tempstring.c_str() << "\t";
        species2.push_back(tempstring);

        getline(infile, tempstring, '\t');
        //std::cout << tempstring.c_str() << "\t";
        etime.push_back(atof(tempstring.c_str()));

        getline(infile, tempstring, '\t');
        //std::cout << tempstring.c_str() << "\n";
        lam_par1.push_back(atof(tempstring.c_str()));

        getline(infile, tempstring, '\t');
        lam_par2.push_back(atof(tempstring.c_str()));

        getline(infile, tempstring, '\t');
        mu_par1.push_back(atof(tempstring.c_str()));

        getline(infile, tempstring, '\n');
        mu_par2.push_back(atof(tempstring.c_str()));

        if (infile.peek() == EOF)
            break;

    }

    infile.close();

    std::cout << "Read a total of " << species1.size() << " events" << std::endl;
    for (std::vector<std::string>::size_type i = 0; i < species1.size(); i++) {
        //std::cout << std::endl << "MRCA of : " <<  species1[i] << "\t" << species2[i] << std::endl;
        if ((species2[i] != "NA") && (species1[i] != "NA")) {

            Node* x = treePtr->getNodeMRCA(species1[i].c_str(), species2[i].c_str());
            if (x  == treePtr->getRoot()) {
                rootEvent->setLamInit(lam_par1[i]);
                rootEvent->setLamShift(lam_par2[i]);
                rootEvent->setMuInit(mu_par1[i]);
                rootEvent->setMuShift(mu_par2[i]);
            } else {
                double deltaT = x->getTime() - etime[i];

                double newmaptime = x->getMapStart() + deltaT;
                //std::cout << std::endl << x->getTime() << "\t" << x->getAnc()->getTime() << std::endl;
                //std::cout << "maptimes: " << x->getMapStart() << "\t" << x->getMapEnd() << "\tcurmap: " << newmaptime << std::endl;
                //std::cout << etime[i] << "\t" << treePtr->getAbsoluteTimeFromMapTime(newmaptime) << std::endl;
                BranchEvent* newEvent =  new BranchEvent(lam_par1[i], lam_par2[i], mu_par1[i],
                        mu_par2[i], x, treePtr, ran, newmaptime);
                newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(newEvent);

                eventCollection.insert(newEvent);
                forwardSetBranchHistories(newEvent);
                treePtr->setMeanBranchSpeciation();
                treePtr->setMeanBranchExtinction();
                //std::cout << newEvent->getAbsoluteTime() << "\t" << newEvent->getLamInit() << "\t" << newEvent->getLamShift() << "\t" << newEvent->getMuInit() << std::endl;
            }

        } else if ((species2[i] == "NA") && (species1[i] != "NA")) {
            Node* x = treePtr->getNodeByName(species1[i].c_str());
            double deltaT = x->getTime() - etime[i];
            double newmaptime = x->getMapStart() + deltaT;
            BranchEvent* newEvent =  new BranchEvent(lam_par1[i], lam_par2[i], mu_par1[i],
                    mu_par2[i], x, treePtr, ran, newmaptime);
            newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(newEvent);

            eventCollection.insert(newEvent);
            forwardSetBranchHistories(newEvent);
            treePtr->setMeanBranchSpeciation();
            treePtr->setMeanBranchExtinction();
            //std::cout << newEvent->getAbsoluteTime() << "\t" << newEvent->getLamInit() << "\t" << newEvent->getLamShift() << "\t" << newEvent->getMuInit() << std::endl;

        } else {
            std::cout << "Error in Model::initializeModelFromEventDataFile" << std::endl;
            exit(1);
        }

        //std::cout << i << "\t" << computeLikelihoodBranches() << "\t" << computeLogPrior() << std::endl;
    }

    std::cout << "Added " << eventCollection.size() <<
         " pre-defined events to tree, plus root event" << std::endl;

    //printEvents();



}



/*
    Adds event to tree based on reference map value
    -adds to branch history set
    -inserts into Model::eventCollection



 */

void Model::addEventToTree(double x)
{

#ifdef ADAPTIVE_MCMC_PROPOSAL

    // For now, the rates of speciation and extinction are set to whatever they should be based
    // on the ancestralNodeEvent
    Node* xnode = treePtr->mapEventToTree(x);
    double atime = treePtr->getAbsoluteTimeFromMapTime(x);
    BranchHistory* bh = xnode->getBranchHistory();
    BranchEvent* be = bh->getAncestralNodeEvent();

    double elapsed = atime - be->getAbsoluteTime();
    double curLam = be->getLamInit() * exp( elapsed * be->getLamShift() );
    double curMu = be->getMuInit() * exp(elapsed * be->getMuShift());
    double curLamShift = be->getLamShift();
    double curMuShift = be->getMuShift();

    // Current parameters.
    // New parameters will be drawn from distributions centered on these.
    // for starters, just draw from exponentials centered on current values for lam and mu
    //  and from a normal for the lamShift parameter.
    double newLam = ran->exponentialRv((double)(1 / curLam));
    double newMu = ran->exponentialRv((double)(1 / curMu));
    double newLambdaShift = ran->normalRv(curLamShift, JUMP_VARIANCE_NORMAL);

    double newMuShift = 0.0;

    _logQratioJump = 0.0; // Set to zero to clear previous values...
    _logQratioJump = ran->lnExponentialPdf((1 / curLam), newLam);
    _logQratioJump += ran->lnExponentialPdf((1 / newMu), newMu);
    _logQratioJump += ran->lnNormalPdf((double)(curLamShift), JUMP_VARIANCE_NORMAL,
                                       newLambdaShift);
    _logQratioJump += ran->lnNormalPdf((double)(curMuShift), JUMP_VARIANCE_NORMAL,
                                       newMuShift);
#else

    // New way:
    // Sample lambda, lambdaInit, mu, muInit parameters from appropriate prior distributions

    double newLam = ran->exponentialRv(_lambdaInitPrior) ;
    double newLambdaShift =  ran->normalRv(0.0, _lambdaShiftPrior);
    double newMu = ran->exponentialRv(_muInitPrior);

    double newMuShift = 0.0;

    // This line would randomly generate a new muShift parameter
    //      but for now we set this to zero
    //double newMuShift = ran->normalRv(0.0, _muShiftPrior);

    // This block computes the jump density for the
    //  addition of new parameters.

    _logQratioJump = 0.0; // Set to zero to clear previous values...
    _logQratioJump = ran->lnExponentialPdf(_lambdaInitPrior, newLam);
    _logQratioJump += ran->lnExponentialPdf(_muInitPrior, newMu);
    _logQratioJump += ran->lnNormalPdf((double)0.0, _lambdaShiftPrior,
                                       newLambdaShift);
    _logQratioJump += ran->lnNormalPdf((double)0.0, _muShiftPrior, newMuShift);

#endif

    // End calculations:: now create event

    BranchEvent* newEvent = new BranchEvent(newLam, newLambdaShift, newMu,
                                            newMuShift, treePtr->mapEventToTree(x), treePtr, ran, x);

    // add the event to the branch history.
    //  ALWAYS done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

    // Addition June17 2012
    lastEventModified = newEvent;

}


/*
 Adds event to tree based on uniform RV
 -adds to branch history set
 -inserts into Model::eventCollection



 */


void Model::addEventToTree(void)
{

    double aa = treePtr->getRoot()->getMapStart();
    double bb = treePtr->getTotalMapLength();
    double x = ran->uniformRv(aa, bb);

#ifdef ADAPTIVE_MCMC_PROPOSAL

    // For now, the rates of speciation and extinction are set to whatever they should be based
    // on the ancestralNodeEvent
    Node* xnode = treePtr->mapEventToTree(x);
    double atime = treePtr->getAbsoluteTimeFromMapTime(x);
    BranchHistory* bh = xnode->getBranchHistory();
    BranchEvent* be = bh->getAncestralNodeEvent();

    double elapsed = atime - be->getAbsoluteTime();
    double curLam = be->getLamInit() * exp( elapsed * be->getLamShift() );
    double curMu = be->getMuInit() * exp(elapsed * be->getMuShift());
    double curLamShift = be->getLamShift();
    double curMuShift = be->getMuShift();

    // Current parameters.
    // New parameters will be drawn from distributions centered on these.
    // for starters, just draw from exponentials centered on current values for lam and mu
    //  and from a normal for the lamShift parameter.
    double newLam = ran->exponentialRv((double)(1 / curLam));
    double newMu = ran->exponentialRv((double)(1 / curMu));
    double newLambdaShift = ran->normalRv(curLamShift, JUMP_VARIANCE_NORMAL);

    double newMuShift = 0.0;

    _logQratioJump = 0.0; // Set to zero to clear previous values...
    _logQratioJump = ran->lnExponentialPdf((1 / curLam), newLam);
    _logQratioJump += ran->lnExponentialPdf((1 / newMu), newMu);
    _logQratioJump += ran->lnNormalPdf((double)(curLamShift), JUMP_VARIANCE_NORMAL,
                                       newLambdaShift);
    _logQratioJump += ran->lnNormalPdf((double)(curMuShift), JUMP_VARIANCE_NORMAL,
                                       newMuShift);

#else

    // New way:
    // Sample lambda, lambdaInit, mu, muInit parameters from appropriate prior distributions

    double newLam = ran->exponentialRv(_lambdaInitPrior) ;
    double newLambdaShift =  ran->normalRv(0.0, _lambdaShiftPrior);
    double newMu = ran->exponentialRv(_muInitPrior);

    double newMuShift = 0.0;

    // This line would randomly generate a new muShift parameter
    //      but for now we set this to zero
    //double newMuShift = ran->normalRv(0.0, _muShiftPrior);

    // This block computes the jump density for the
    //  addition of new parameters.

    _logQratioJump = 0.0; // Set to zero to clear previous values...
    _logQratioJump = ran->lnExponentialPdf(_lambdaInitPrior, newLam);
    _logQratioJump += ran->lnExponentialPdf(_muInitPrior, newMu);
    _logQratioJump += ran->lnNormalPdf((double)0.0, _lambdaShiftPrior,
                                       newLambdaShift);
    _logQratioJump += ran->lnNormalPdf((double)0.0, _muShiftPrior, newMuShift);

#endif

    BranchEvent* newEvent = new BranchEvent(newLam, newLambdaShift, newMu,
                                            newMuShift, treePtr->mapEventToTree(x), treePtr, ran, x);



    newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);
    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

    // Addition June17 2012
    lastEventModified = newEvent;

}

void Model::printEvents(void)
{

    // for each event:
    //  print:  maptime
    //          nodeptr
    //
    int n_events = (int)eventCollection.size();
    std::cout << "N_events: " << n_events << std::endl;
    int counter = 1;
    for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); i++)
        std::cout << "event " << counter++ << "\tAddress: " << (*i) << "\t" <<
             (*i)->getMapTime() << "\tNode: " << (*i)->getEventNode() << std::endl << std::endl;


}

BranchEvent* Model::chooseEventAtRandom(void)
{

    int n_events = (int)eventCollection.size();
    if (n_events == 0) {
        return NULL;
        //should ultimately throw exception here.

    } else {
        int ctr = 0;
        double xx = ran->uniformRv();
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

void Model::eventLocalMove(void)
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
		double step = ran->uniformRv(0, _scale) - 0.5*_scale;
		
        chosenEvent->moveEventLocal(step); // move event
        chosenEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
            chosenEvent);

        // Get last event from the theEventNode
        //      forward set its history
        //  Then go to the "moved" event and forward set its history


        forwardSetBranchHistories(previousEvent);
        forwardSetBranchHistories(chosenEvent);

    }
    // else no events to move
    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

}

void Model::eventGlobalMove(void)
{

    if (getNumberOfEvents() > 0) {
        BranchEvent* chosenEvent = chooseEventAtRandom();

        // this is the event preceding the chosen event: histories should be set forward from here..
        BranchEvent* previousEvent =
            chosenEvent->getEventNode()->getBranchHistory()->getLastEvent(chosenEvent);

        //Node* theEventNode = chosenEvent->getEventNode();

        // private variable
        lastEventModified = chosenEvent;

        chosenEvent->getEventNode()->getBranchHistory()->popEventOffBranchHistory(
            chosenEvent);
        chosenEvent->moveEventGlobal(); // move event
        chosenEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(
            chosenEvent);

        // Get previous event from the theEventNode
        //      forward set its history
        //  Then go to the "moved" event and forward set its history

        forwardSetBranchHistories(previousEvent);
        forwardSetBranchHistories(chosenEvent);

    }
    //std::cout << "leave globalMove" << std::endl;

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

}

// used to reset position of event if move is rejected

void Model::revertMovedEventToPrevious(void)
{

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
    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

}



// Recursively count the number of events in the branch histories
int Model::countEventsInBranchHistory(Node* p)
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

void Model::deleteEventFromTree(BranchEvent* be)
{

    if (be == rootEvent) {
        std::cout << "Can't delete root event" << std::endl;
        exit(1);
    } else {
        // erase from branch history:
        Node* currNode = (be)->getEventNode();

        //get event downstream of i
        BranchEvent* newLastEvent = currNode->getBranchHistory()->getLastEvent(be);

        lastDeletedEventMapTime = (be)->getMapTime();
        _lastDeletedEventLambdaInit = (be)->getLamInit();
        _lastDeletedEventLambdaShift = (be)->getLamShift();

        _lastDeletedEventMuInit = (be)->getMuInit();
        _lastDeletedEventMuShift = (be)->getMuShift();


        _logQratioJump = 0.0; // Set to zero to clear previous values...
        _logQratioJump = ran->lnExponentialPdf(_lambdaInitPrior,
                                               _lastDeletedEventLambdaInit);
        _logQratioJump += ran->lnExponentialPdf(_muInitPrior, _lastDeletedEventMuInit);
        _logQratioJump += ran->lnNormalPdf((double)0.0, _lambdaShiftPrior,
                                           _lastDeletedEventLambdaShift);
        _logQratioJump += ran->lnNormalPdf((double)0.0, _muShiftPrior,
                                           _lastDeletedEventMuShift);


        currNode->getBranchHistory()->popEventOffBranchHistory((be));

        eventCollection.erase(be);

        // delete from global node set
        delete (be);
        //std::cout << "deleted..." << std::endl;

        forwardSetBranchHistories(newLastEvent);
    }


    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

}




void Model::deleteRandomEventFromTree(void)
{

    //std::cout << std::endl << std::endl << "START Delete: " << std::endl;
    //printBranchHistories(treePtr->getRoot());

    // can only delete event if more than root node present.
    int n_events = (int)eventCollection.size();

    if (eventCollection.size() > 0) {
        int counter = 0;
        double xx = ran->uniformRv();
        int chosen = (int)(xx * (double)n_events);

        for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
                i != eventCollection.end(); i++) {
            if (counter++ == chosen) {

                // erase from branch history:
                Node* currNode = (*i)->getEventNode();

                //get event downstream of i
                BranchEvent* newLastEvent = currNode->getBranchHistory()->getLastEvent((*i));

                lastDeletedEventMapTime = (*i)->getMapTime();
                //lastDeletedEventBeta = (*i)->getBeta();

                _lastDeletedEventLambdaInit = (*i)->getLamInit();
                _lastDeletedEventLambdaShift = (*i)->getLamShift();

                _lastDeletedEventMuInit = (*i)->getMuInit();
                _lastDeletedEventMuShift = (*i)->getMuShift();

#ifdef ADAPTIVE_MCMC_PROPOSAL


                double atime = treePtr->getAbsoluteTimeFromMapTime(lastDeletedEventMapTime);

                double elapsed = atime - newLastEvent->getAbsoluteTime();
                double newLam = newLastEvent->getLamInit() * exp( elapsed *
                                newLastEvent->getLamShift() );
                double newMu = newLastEvent->getMuInit() * exp(elapsed *
                               newLastEvent->getMuShift());
                double newLamShift = newLastEvent->getLamShift();
                double newMuShift = newLastEvent->getMuShift();

                _logQratioJump = 0.0; // Set to zero to clear previous values...
                _logQratioJump = ran->lnExponentialPdf((1 / newLam),
                                                       _lastDeletedEventLambdaInit);
                _logQratioJump += ran->lnExponentialPdf((1 / newMu), _lastDeletedEventMuInit);
                _logQratioJump += ran->lnNormalPdf((double)(newLamShift), JUMP_VARIANCE_NORMAL,
                                                   _lastDeletedEventLambdaShift);
                _logQratioJump += ran->lnNormalPdf((double)(newMuShift), JUMP_VARIANCE_NORMAL,
                                                   _lastDeletedEventMuShift);



#else

                // This block computes the jump density for the
                //  deletion of new parameters.


                _logQratioJump = 0.0; // Set to zero to clear previous values...
                _logQratioJump = ran->lnExponentialPdf(_lambdaInitPrior,
                                                       _lastDeletedEventLambdaInit);
                _logQratioJump += ran->lnExponentialPdf(_muInitPrior, _lastDeletedEventMuInit);
                _logQratioJump += ran->lnNormalPdf((double)0.0, _lambdaShiftPrior,
                                                   _lastDeletedEventLambdaShift);
                _logQratioJump += ran->lnNormalPdf((double)0.0, _muShiftPrior,
                                                   _lastDeletedEventMuShift);

#endif

                //std::cout << (*i) << std::endl;

                currNode->getBranchHistory()->popEventOffBranchHistory((*i));

                eventCollection.erase(i);

                // delete from global node set
                delete (*i);

                forwardSetBranchHistories(newLastEvent);

                break;
            }
        }
    }

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

}


void Model::restoreLastDeletedEvent(void)
{


    // Use constructor for speciation and extinction

    BranchEvent* newEvent = new BranchEvent((double)0.0, (double)0.0, (double)0.0,
                                            (double)0.0, treePtr->mapEventToTree(lastDeletedEventMapTime), treePtr, ran,
                                            lastDeletedEventMapTime);

    newEvent->setLamInit(_lastDeletedEventLambdaInit);
    newEvent->setLamShift(_lastDeletedEventLambdaShift);
    newEvent->setMuInit(_lastDeletedEventMuInit);
    newEvent->setMuShift(_lastDeletedEventMuShift);



    // add the event to the branch history.
    //  ALWAYS done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(newEvent);

    eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

}



void Model::changeNumberOfEventsMH(void)
{


    // Get old prior density of the data:
    double oldLogPrior = computeLogPrior();
    double newLogPrior = 0.0;
    int currState = (int)eventCollection.size();
    int proposedState = 0;
    bool acceptMove = false;

    // Propose gains & losses equally if not on boundary (n = 0) events:

    // Current number of events on the tree, not counting root state:
    double K = (double)(eventCollection.size());

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
        double likBranches = 0.0;
        double PropLnLik = likBranches;

#else

        addEventToTree();

        treePtr->setMeanBranchSpeciation();
        treePtr->setMeanBranchExtinction();

        double likBranches = computeLikelihoodBranches();
        double PropLnLik = likBranches;



#endif
        // Prior density on all parameters
        newLogPrior = computeLogPrior();


        // Prior ratio is eventRate / (k + 1)
        // but now, eventCollection.size() == k + 1
        //  because event has already been added.
        // Here HR is just the prior ratio

        double logHR = log(eventLambda) - log(K + 1.0);

        // Now add log qratio

        logHR += log(qratio);

        double likeRatio = PropLnLik - getCurrLnLBranches();

        logHR += likeRatio;

        // Now for prior:
        logHR += (newLogPrior - oldLogPrior);

        // Now for jumping density of the bijection between parameter spaces:
        logHR -= _logQratioJump;

        if (std::isinf(likeRatio) ) {

        } else
            acceptMove = acceptMetropolisHastings(logHR);




        if (acceptMove) {
            //std::cout << "gaining event in changeNumberOfEventsMH " << std::endl;
            //std::cout << "gainaccept" << computeLikelihoodBranches()  << std::endl;

            //addEventToTree();
            //std::cout << "Calliing isValid from ChangeNumberEvents::gain" << std::endl;

            bool isValidConfig = isEventConfigurationValid(lastEventModified);

            if (isValidConfig) {
                // Update accept/reject statistics
                acceptCount++;
                acceptLast = 1;
            } else {
                // Need to get rid of event that was just gained...
                //std::cout << "Invalid event config from addEventToTree - deleting." << std::endl;
                deleteEventFromTree(lastEventModified);
                rejectCount++;
                acceptLast = 0;

            }


        } else {
            //std::cout << "gainreject" << computeLikelihoodBranches() << std::endl;

            // Delete event.
            deleteEventFromTree(lastEventModified);

            rejectCount++;
            acceptLast = 0;
        }





    } else {
        //std::cout << "loss: initial LH: " << computeLikelihoodBranches() << std::endl;
        deleteRandomEventFromTree();
        //std::cout << "loss: LH after deleteRandomEvent" << computeLikelihoodBranches() << std::endl;

        proposedState = currState - 1;

#ifdef NO_DATA
        double likBranches = 0.0;
        double PropLnLik = likBranches;

#else

        treePtr->setMeanBranchSpeciation();
        treePtr->setMeanBranchExtinction();

        double likBranches = computeLikelihoodBranches();
        double PropLnLik = likBranches;
        double likTraits = 0.0;

#endif

        // Prior density on all parameters
        newLogPrior = computeLogPrior();

        double qratio = 1.0; // if loss, can only be qratio of 1.0

        if (K  == 1)
            qratio = 2.0;

        // This is probability of going from k to k-1
        // So, prior ratio is (k / eventRate)

        // First get prior ratio:
        double logHR = log(K) - log(eventLambda);

        // Now correct for proposal ratio:
        logHR += log(qratio);

        double likeRatio = PropLnLik - getCurrLnLBranches();


        logHR += likeRatio;

        // Now for prior:
        logHR += (newLogPrior - oldLogPrior);

        // Now for jumping density of the bijection between parameter spaces:
        logHR += _logQratioJump;

        if (std::isinf(likeRatio) ) {

        } else
            acceptMove = acceptMetropolisHastings(logHR);


        //std::cout << "loss: " << acceptMove << "\t" << PropLnLik << "\tLT " << getCurrLnLTraits() + getCurrLnLBranches() << std::endl;

        if (acceptMove) {
            //std::cout << "loss accept, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;
            setCurrLnLTraits(likTraits);

            setCurrLnLBranches(likBranches);

            acceptCount++;
            acceptLast = 1;


        } else {
            //std::cout << "loss reject, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;

            restoreLastDeletedEvent();

            // speciation-extinction rates on branches automatically updated after restoreLastDeletedEvent()

            //std::cout << "loss reject restored, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;
            rejectCount++;
            acceptLast = 0;
        }

    }

    //std::cout << currState << "\t" << proposedState << "\tAcc: " << acceptMove << "\tOP: " << oldLogPrior;
    //std::cout << "\tNP: " << newLogPrior << "\tqratio: " << _logQratioJump << std::endl;

    incrementGeneration();

}

void Model::moveEventMH(void)
{


    if (eventCollection.size() > 0) {

        double localMoveProb = _localGlobalMoveRatio / (1 + _localGlobalMoveRatio);

        bool isLocalMove = (ran->uniformRv() <= localMoveProb);
        //std::cout << "is local: " << isLocalMove << std::endl;

        if (isLocalMove) {
            // Local move, with event drawn at random
            eventLocalMove();

#ifdef NO_DATA
            double likBranches = 0;
            double PropLnLik = likBranches;
#else
            treePtr->setMeanBranchSpeciation();
            treePtr->setMeanBranchExtinction();

            double likBranches = computeLikelihoodBranches();
            double PropLnLik = likBranches;



#endif

            double likeRatio = PropLnLik - getCurrLnLBranches();
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

                setCurrLnLBranches(likBranches);
                acceptCount++;
                acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                treePtr->setMeanBranchSpeciation();
                treePtr->setMeanBranchExtinction();

                rejectCount++;
                acceptLast = 0;
            }



        } else {
            // global move, event drawn at random
            eventGlobalMove();
            //std::cout << "successful global move" << std::endl;
#ifdef NO_DATA
            double likBranches = 0;
            double PropLnLik = likBranches;
#else

            treePtr->setMeanBranchSpeciation();
            treePtr->setMeanBranchExtinction();

            double likBranches = computeLikelihoodBranches();
            double PropLnLik = likBranches;

#endif

            double likeRatio = PropLnLik - getCurrLnLBranches();
            double logHR = likeRatio;


            bool acceptMove = false;
            bool isValid = false;

            isValid = isEventConfigurationValid(lastEventModified);

            if (std::isinf(likeRatio) ) {

            } else if (isValid)
                acceptMove = acceptMetropolisHastings(logHR);
            else {

            }

            if (acceptMove == true) {

                setCurrLnLBranches(likBranches);
                acceptCount++;
                acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                treePtr->setMeanBranchSpeciation();
                treePtr->setMeanBranchExtinction();

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




void Model::updateLambdaInitMH(void)
{

    //int n_events = eventCollection.size() + 1;
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    BranchEvent* be = rootEvent;

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = (*myIt);
    }

    double oldRate = be->getLamInit();
    double cterm = exp( _updateLambdaInitScale * (ran->uniformRv() - 0.5) );

    be->setLamInit(cterm * oldRate);

#ifdef NO_DATA
    double PropLnLik = 0;
#else

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

        
    double logPriorRatio = cprior->lambdaInitPrior(be->getLamInit());
    logPriorRatio -= cprior->lambdaInitPrior(oldRate);
    
    
    double LogProposalRatio = log(cterm);

    double likeRatio = PropLnLik - getCurrLnLBranches();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrLnLBranches(PropLnLik);
        acceptCount++;
        acceptLast = 1;
    } else {

        // revert to previous state
        be->setLamInit(oldRate);

        treePtr->setMeanBranchSpeciation();
        treePtr->setMeanBranchExtinction();

        acceptLast = 0;
        rejectCount++;
    }

    incrementGeneration();


}

void Model::updateLambdaShiftMH(void)
{

    //int n_events = eventCollection.size() + 1;
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    BranchEvent* be = rootEvent;

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = (*myIt);
    }

    double oldLambdaShift = be->getLamShift();
    double newLambdaShift = oldLambdaShift + ran->normalRv((double)0.0,
                            _updateLambdaShiftScale);
    be->setLamShift(newLambdaShift);


#ifdef NO_DATA
    double PropLnLik = 0;
#else

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

    double  logPriorRatio = cprior->lambdaShiftPrior(newLambdaShift);
    logPriorRatio -= cprior->lambdaShiftPrior(oldLambdaShift);
    
/*
    double  logPriorRatio = ran->lnNormalPdf((double)0.0,
                            sttings->getLambdaShiftPrior(), newLambdaShift);
    logPriorRatio -= ran->lnNormalPdf((double)0.0, sttings->getLambdaShiftPrior(),
                                      oldLambdaShift);
*/
    
    double LogProposalRatio = 0.0;

    double likeRatio = PropLnLik - getCurrLnLBranches();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);



    if (acceptMove == true) {

        setCurrLnLBranches(PropLnLik);
        acceptCount++;
        acceptLast = 1;
    } else {

        // revert to previous state
        be->setLamShift(oldLambdaShift);

        treePtr->setMeanBranchSpeciation();
        treePtr->setMeanBranchExtinction();

        acceptLast = 0;
        rejectCount++;
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
void Model::updateTimeVariablePartitionsMH(void)
{

    //int n_events = eventCollection.size() + 1;
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    BranchEvent* be = rootEvent;

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


void Model::updateMuInitMH(void)
{

    //int n_events = eventCollection.size() + 1;
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    BranchEvent* be = rootEvent;

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = (*myIt);
    }

    double oldRate = be->getMuInit();
    double cterm = exp( _updateMuInitScale * (ran->uniformRv() - 0.5) );

    be->setMuInit(cterm * oldRate);

#ifdef NO_DATA
    double PropLnLik = 0;
#else

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

    double logPriorRatio = cprior->muInitPrior(be->getMuInit());
    logPriorRatio -= cprior->muInitPrior(oldRate);

    double LogProposalRatio = log(cterm);

    double likeRatio = PropLnLik - getCurrLnLBranches();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrLnLBranches(PropLnLik);
        acceptCount++;
        acceptLast = 1;
    } else {

        // revert to previous state
        be->setMuInit(oldRate);

        treePtr->setMeanBranchSpeciation();
        treePtr->setMeanBranchExtinction();

        acceptLast = 0;
        rejectCount++;
    }

    incrementGeneration();


}


void Model::updateMuShiftMH(void)
{

    //int n_events = eventCollection.size() + 1;
    int toUpdate = ran->sampleInteger(0, (int)eventCollection.size());
    BranchEvent* be = rootEvent;

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = (*myIt);
    }

    double oldMuShift = be->getMuShift();
    double newMuShift = oldMuShift + ran->normalRv((double)0.0,
                        _updateMuShiftScale);

    be->setMuShift(newMuShift);

#ifdef NO_DATA
    double PropLnLik = 0;
#else

    treePtr->setMeanBranchSpeciation();
    treePtr->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

    
    double logPriorRatio = cprior->muShiftPrior(newMuShift);
    logPriorRatio -= cprior->muShiftPrior(oldMuShift);
    
    

    double LogProposalRatio = 0.0;

    double likeRatio = PropLnLik - getCurrLnLBranches();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;



    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrLnLBranches(PropLnLik);
        acceptCount++;
        acceptLast = 1;

    } else {

        // revert to previous state
        be->setMuShift(oldMuShift);

        treePtr->setMeanBranchSpeciation();
        treePtr->setMeanBranchExtinction();

        acceptLast = 0;
        rejectCount++;
    }

    incrementGeneration();


}




/*

 Metropolis-Hastings step to update Poisson event rate.
 Note that changing this rate does not affect the likelihood,
 so the priors and qratio determine acceptance rate.

 */
void Model::updateEventRateMH(void)
{


    double oldEventRate = getEventRate();

    double cterm = exp( _updateEventRateScale * (ran->uniformRv() - 0.5) );
    setEventRate(cterm * oldEventRate);

    
    double logPriorRatio = cprior->poissonRatePrior(getEventRate());
    logPriorRatio -= cprior->poissonRatePrior(oldEventRate);
    

    double logProposalRatio = log(cterm);


    // Experimental code:
    // Sample new event rate from prior directly with each step.
    //double newEventRate = ran->exponentialRv(_poissonRatePrior);
    //setEventRate(newEventRate);
    //double LogPriorRatio = 0.0;
    //double logProposalRatio = 1.0;

    double logHR = logPriorRatio + logProposalRatio;
    const bool acceptMove = acceptMetropolisHastings(logHR);


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

}


double Model::computeLikelihoodBranches(void)
{

    return computeLikelihoodBranchesByInterval();
}


double Model::computeLikelihoodBranchesByInterval(void)
{


    double LnL = 0.0;


    if (sttings->getSampleFromPriorOnly())
        return 0.0;

    int numNodes = treePtr->getNumberOfNodes();

    // LEft and right extinction probabilities for root node
    double rootEleft = 0.0;
    double rootEright = 0.0;


    for (int i = 0; i < numNodes; i++) {
        Node* xnode = treePtr->getNodeFromDownpassSeq(i);
        if (xnode->getLfDesc() != NULL && xnode->getRtDesc() != NULL) {
            // NOT tip, but MUST ultimately be root.

            // Do left descendant:
            Node* ldesc = xnode->getLfDesc();

            double lDinit = ldesc->getDinit();
            double lEinit = ldesc->getEinit();
            double starttime = ldesc->getBrlen();
            double endtime = ldesc->getBrlen();

            double LtotalT = 0.0; // what is this doing?

            double meanRateAtBranchBase = 0.0;
            double curLam = 0.0;

            while (starttime > 0) {
                //std::cout << starttime << "\t" << endtime << std::endl;
                starttime -= _segLength;
                if (starttime < 0)
                    starttime = 0.0;
                double deltaT = endtime - starttime;

                LtotalT += deltaT;

                curLam = ldesc->computeSpeciationRateIntervalRelativeTime(starttime, endtime);

                double curMu = ldesc->computeExtinctionRateIntervalRelativeTime(starttime,
                               endtime);

                double numL = 0.0;
                double denomL = 0.0;


                numL = (exp( deltaT * (curMu - curLam)) * lDinit * ((curLam - curMu) *
                        (curLam - curMu) ) );
                denomL = ( curLam - (lEinit * curLam) + (exp(deltaT * (curMu - curLam)) *
                           (lEinit * curLam - curMu)));

                lDinit = (numL / (denomL * denomL));
                LnL += log(lDinit);
                lDinit = 1.0;


                double Enum = (1 - lEinit) * (curLam - curMu);
                double Edenom =  (1 - lEinit) * curLam - (exp((curMu - curLam) * (deltaT))) *
                                 (curMu - curLam * lEinit);

                lEinit = 1.0 - (Enum / Edenom);


                endtime = starttime; // reset starttime to old endtime
            }

            // this is to get node speciation rate using approximations
            //   to correspond to fact that branch likelihoods themselves are computed
            //      using approximations.

            meanRateAtBranchBase = curLam / 2;
            curLam = 0.0;

            // Setting extinction prob at root node IF xnode is the root
            if (xnode == treePtr->getRoot())
                rootEleft = lEinit;

            // Compute speciation for right descendant
            // Do right descendant:
            Node* rdesc = xnode->getRtDesc();

            double rDinit = rdesc->getDinit();
            double rEinit = rdesc->getEinit();

            starttime = rdesc->getBrlen();
            endtime = rdesc->getBrlen();

            double RtotalT = 0.0;

            while (starttime > 0) {

                starttime -= _segLength;
                if (starttime < 0)
                    starttime = 0.0;
                double deltaT = endtime - starttime;

                RtotalT += deltaT;

                curLam = rdesc->computeSpeciationRateIntervalRelativeTime(starttime, endtime);

                double curMu = rdesc->computeExtinctionRateIntervalRelativeTime(starttime,
                               endtime);

                double numL = 0.0;
                double denomL = 0.0;

                numL = (exp( deltaT * (curMu - curLam)) * rDinit * ((curLam - curMu) *
                        (curLam - curMu) ) );
                denomL = ( curLam - (rEinit * curLam) + (exp(deltaT * (curMu - curLam)) *
                           (rEinit * curLam - curMu)));

                rDinit = (numL / (denomL * denomL));
                LnL += log(rDinit);
                rDinit = 1.0;

                double Enum = 0.0;
                double Edenom = 0.0;

                Enum = (1 - rEinit) * (curLam - curMu);
                Edenom =  (1 - rEinit) * curLam - (exp((curMu - curLam) * (deltaT))) *
                          (curMu - curLam * rEinit);


                rEinit = 1.0 - (Enum / Edenom);



                endtime = starttime; // reset starttime to old endtime


            }

            meanRateAtBranchBase += curLam / 2;

            // ########################## What to use as Einit for start of NEXT downstream branch?
            // At this point, lEinit is actually the lEinit for the parent node:
            //  Here, we will  (as of 9.1.2012) arbitrarily take this to be the LEFT extinction value:
            xnode->setEinit(lEinit);

            // ######## But alternatively, we could do:
            // Like the above, but now we randomly resolve this. We choose at RANDOM whether to use the right or left Evalue.
            // as we don't know which descendant represents the "parent" state.

            /*          ***************
            if (ran->uniformRv() <= 0.5){
                xnode->setEinit(lEinit);
            }else{
                xnode->setEinit(rEinit);
            }
                        ****************        */


            // Clearly a problem if extinction values approaching/equaling 1.0
            // If so, set to -Inf, leading to automatic rejection of state

            if ((lEinit > MAX_E_PROB) || (rEinit > MAX_E_PROB)) {
                //std::cout << xnode << "\t" << lEinit << "\t" << rEinit << std::endl;
                return -INFINITY;
            }


            if (xnode == treePtr->getRoot())
                rootEright = rEinit;
            // rDinit at this point should be FINAL value:
            // save as new variable, to keep clear:

            /* SHould be abele to ignore all of these calculations for the root node:
             Must also compute speciation rate for focal node. THis is a critical step.

             Since we are using approximations for the calculations on branches, we should set node speciation
             rate to be equivalent. Currently, I am not doing this - just computing exact rates
             at nodes.
            */

            if (xnode != treePtr->getRoot()) {

                // Does not include root node, so it is conditioned on basal speciation event occurring:

                LnL  += log(xnode->getNodeLambda());

                //LnL += log(meanRateAtBranchBase);



                xnode->setDinit(1.0);


            }



        } // IF not tip

    } // FOR each node in set



    // 09.15.2012
    // To CONDITION, uncomment the line below:
    // Or, if UNCOMMENTED, comment the line to NOT condition on survival
    LnL -= (log(1 - rootEleft) + log(1 -
                                     rootEright)); // replacement to above for condiioning.


    return LnL;
}



/* Only works on speciation + extinction
 */

double Model::computeLogPrior(void)
{

    double logPrior = 0.0;

    logPrior += cprior->lambdaInitPrior(rootEvent->getLamInit());
    logPrior += cprior->lambdaShiftPrior(rootEvent->getLamShift());
    logPrior += cprior->muInitPrior(rootEvent->getMuInit());
    logPrior += cprior->muShiftPrior(rootEvent->getMuShift());
    
    int ctr = 0;

    for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); i++) {

        logPrior += cprior->lambdaInitPrior((*i)->getLamInit());
        logPrior += cprior->lambdaShiftPrior((*i)->getLamShift());
        logPrior += cprior->muInitPrior((*i)->getMuInit());
        logPrior += cprior->muShiftPrior((*i)->getMuShift());
        
        ctr++;

    }

    // Here's prior density on the event rate:
    
    logPrior += cprior->poissonRatePrior(getEventRate());
    
    return logPrior;

}





bool Model::acceptMetropolisHastings(const double lnR)
{
    const double r = safeExponentiation(Model::mhColdness * lnR);
    return (ran->uniformRv() < r);
}




void Model::initializeBranchHistories(Node* x)
{
    //std::cout << x << std::endl;
    x->getBranchHistory()->setNodeEvent(rootEvent);

    if (x->getAnc() != NULL)
        x->getBranchHistory()->setAncestralNodeEvent(rootEvent);

    if (x->getLfDesc() != NULL)
        initializeBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        initializeBranchHistories(x->getRtDesc());


}



void Model::printStartAndEndEventStatesForBranch(Node* x)
{

    if (x != treePtr->getRoot())
        std::cout << "Node: " << x << "\tAnc: " <<
             x->getBranchHistory()->getAncestralNodeEvent() << "\tevent: " <<
             x->getBranchHistory()->getNodeEvent() << std::endl;

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

void Model::forwardSetBranchHistories(BranchEvent* x)
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


void Model::forwardSetHistoriesRecursive(Node* p)
{

    // Get event that characterizes parent node
    BranchEvent* lastEvent = p->getAnc()->getBranchHistory()->getNodeEvent();
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




void Model::printBranchHistories(Node* x)
{

    if (x != treePtr->getRoot()) {
        std::cout << "Node: " << x;
        std::cout << "\t#Events: " << x->getBranchHistory()->getNumberOfBranchEvents() <<
             "\tStart: ";
        std::cout << x->getBranchHistory()->getAncestralNodeEvent() << "\tEnd: ";
        std::cout << x->getBranchHistory()->getNodeEvent() << std::endl;

    }
    if (x->getLfDesc() != NULL)
        printBranchHistories(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        printBranchHistories(x->getRtDesc());



}



double  Model::getMHacceptanceRate(void)
{

    double arate = (double)acceptCount / ((double)acceptCount +
                                          (double)rejectCount);

    return arate;

}


void  Model::resetMHacceptanceParameters(void)
{
    acceptCount = 0;
    rejectCount = 0;

}



BranchEvent* Model::getEventByIndex(int x)
{

    //int ctr = 0;
    std::set<BranchEvent*>::iterator myIt = eventCollection.begin();
    for (int i = 0; i <= x; i++)
        myIt++;

    return (*myIt);
}


// adding log contrasts here.


void Model::printExtinctionParams(void)
{

    if (eventCollection.size() > 0) {
        for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
                i != eventCollection.end(); i++)
            std::cout << (*i) << "\t" << (*i)->getMuInit() << "\t" << (*i)->getMuShift() << std::endl;

    }
    std::cout << rootEvent << "\t" << rootEvent->getMuInit() << "\t" <<
         rootEvent->getMuShift() << std::endl << std::endl;


}


/*
 Model::countTimeVaryingRatePartitions

    -counts number of time-varying rate partitions

 */
int Model::countTimeVaryingRatePartitions(void)
{

    int count = 0;
    count += (int)rootEvent->getIsEventTimeVariable();
    for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
            i != eventCollection.end(); i++)
        count += (int)(*i)->getIsEventTimeVariable();
    return count;
}


/*
    Write event data to file for all events "on" tree
        at a given point in the MCMC chain


 */

void Model::getEventDataString(std::stringstream& ss)
{

    ss << getGeneration() << ",";


    BranchEvent* be = rootEvent;
    Node* xl = treePtr->getRoot()->getRandomLeftTipNode();
    Node* xr = treePtr->getRoot()->getRandomRightTipNode();
    ss << xl->getName() << "," << xr->getName() << "," << be->getAbsoluteTime() <<
       ",";
    ss << be->getLamInit() << "," << be->getLamShift() << "," << be->getMuInit() <<
       ",";
    ss << be->getMuShift();




    if (eventCollection.size() > 0) {
        for (std::set<BranchEvent*>::iterator i = eventCollection.begin();
                i != eventCollection.end(); i++) {

            ss << "\n" << getGeneration() << ",";
            be = (*i);
            if (be->getEventNode()->getLfDesc() == NULL)
                ss << be->getEventNode()->getName() << "," << "NA" << ",";

            else {
                Node* xl = be->getEventNode()->getRandomLeftTipNode();
                Node* xr = be->getEventNode()->getRandomRightTipNode();
                ss << xl->getName() << "," << xr->getName() << ",";
            }
            ss << be->getAbsoluteTime() << ",";
            ss << be->getLamInit() << "," << be->getLamShift() << "," << be->getMuInit()  <<
               "," << be->getMuShift();

        }

    }
}


bool Model::isEventConfigurationValid(BranchEvent* be)
{
    //std::cout << "enter isEventConfigValid" << std::endl;
    bool isValidConfig = false;

    if (be->getEventNode() == treePtr->getRoot()) {
        Node* rt = treePtr->getRoot()->getRtDesc();
        Node* lf = treePtr->getRoot()->getLfDesc();
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

        if (anc == treePtr->getRoot())
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




double safeExponentiation(double x)
{

    if (x > 0.0)
        return 1.0;
    else if (x < -100.0)
        return 0.0;
    else
        return exp(x);

}


void Model::debugLHcalculation(void)
{
    std::cout << "This does not currently support anything" << std::endl;

}









