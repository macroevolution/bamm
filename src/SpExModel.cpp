#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>

#include "SpExModel.h"
#include "Model.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "SpExBranchEvent.h"
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


SpExModel::SpExModel(MbRandom* ranptr, Tree* tp, Settings* sp, Prior* pr) :
    Model(ranptr, tp, sp, pr)
{
    // Initialize MCMC proposal/tuning parameters
    _updateLambdaInitScale = _settings->getUpdateLambdaInitScale();
    _updateMuInitScale = _settings->getUpdateMuInitScale();
    _updateLambdaShiftScale = _settings->getUpdateLambdaShiftScale();
    _updateMuShiftScale = _settings->getUpdateMuShiftScale();

    // Initial values
    _lambdaInit0 = _settings->getLambdaInit0();
    _lambdaShift0 = _settings->getLambdaShift0();
    _muInit0 = _settings->getMuInit0();
    _muShift0 = _settings->getMuShift0();

    // For Poisson process
    _lambdaInitPrior = _settings->getLambdaInitPrior();
    _lambdaShiftPrior = _settings->getLambdaShiftPrior();
    _muInitPrior = _settings->getMuInitPrior();
    _muShiftPrior = _settings->getMuShiftPrior();

    // Parameter for splitting branch into pieces for numerical computation
    _segLength = _settings->getSegLength() * _tree->maxRootToTipLength();
   
    BranchEvent* x =  new SpExBranchEvent
        (_lambdaInit0, _lambdaShift0, _muInit0, _muShift0,
            _tree->getRoot(), _tree, _rng, 0);
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the_rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event
    forwardSetBranchHistories(_rootEvent);

    _tree->setMeanBranchSpeciation();
    _tree->setMeanBranchExtinction();

    // Initialize by previous event histories
    if (_settings->getLoadEventData()) {
        log() << "\nLoading model data from file.\n";
        initializeModelFromEventDataFile();
    }

    setCurrLnLBranches(computeLikelihoodBranches());

    log() << "\nInitial log-likelihood: " << getCurrLnLBranches() << "\n";
    if (_settings->getSampleFromPriorOnly())
        log() << "Note that you have chosen to sample from prior only.\n";
}


SpExModel::~SpExModel(void)
{
    for (std::set<BranchEvent*>::iterator it = _eventCollection.begin();
            it != _eventCollection.end(); it++)
        delete (*it);
}


void SpExModel::readModelSpecificParameters(std::ifstream &inputFile)
{
    inputFile >> _readLambdaInit;
    inputFile >> _readLambdaShift;
    inputFile >> _readMuInit;
    inputFile >> _readMuShift;
}


void SpExModel::setRootEventWithReadParameters()
{
    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    rootEvent->setLamInit(_readLambdaInit);
    rootEvent->setLamShift(_readLambdaShift);
    rootEvent->setMuInit(_readMuInit);
    rootEvent->setMuShift(_readMuShift);
}


BranchEvent* SpExModel::newBranchEventWithReadParameters(Node* x, double time)
{
    return new SpExBranchEvent(_readLambdaInit, _readLambdaShift,
        _readMuInit, _readMuShift, x, _tree, _rng, time);
}


void SpExModel::setMeanBranchParameters()
{
    _tree->setMeanBranchSpeciation();
    _tree->setMeanBranchExtinction();
}


BranchEvent* SpExModel::newBranchEventWithRandomParameters(double x)
{
    double newLam = _prior->generateLambdaInitFromPrior();
    double newLambdaShift = _prior->generateLambdaShiftFromPrior();
    double newMu = _prior->generateMuInitFromPrior();
    double newMuShift = _prior->generateMuShiftFromPrior();
 
    // Computes the jump density for the addition of new parameters.
    _logQRatioJump = 0.0;    // Set to zero to clear previous values
    _logQRatioJump = _prior->lambdaInitPrior(newLam);
    _logQRatioJump += _prior->lambdaShiftPrior(newLambdaShift);
    _logQRatioJump += _prior->muInitPrior(newMu);
    _logQRatioJump += _prior->muShiftPrior(newMuShift);
    
    return new SpExBranchEvent(newLam, newLambdaShift, newMu,
        newMuShift, _tree->mapEventToTree(x), _tree, _rng, x);
}


void SpExModel::setDeletedEventParameters(BranchEvent* be)
{
    SpExBranchEvent* event = static_cast<SpExBranchEvent*>(be);

    _lastDeletedEventLambdaInit = event->getLamInit();
    _lastDeletedEventLambdaShift = event->getLamShift();
    _lastDeletedEventMuInit = event->getMuInit();
    _lastDeletedEventMuShift = event->getMuShift();
}


double SpExModel::calculateLogQRatioJump()
{
    double _logQRatioJump = 0.0;
    
    _logQRatioJump = _prior->lambdaInitPrior(_lastDeletedEventLambdaInit);
    _logQRatioJump += _prior->lambdaShiftPrior(_lastDeletedEventLambdaShift);
    _logQRatioJump += _prior->muInitPrior(_lastDeletedEventMuInit);
    _logQRatioJump += _prior->muShiftPrior(_lastDeletedEventMuShift);

    return _logQRatioJump;
}


void SpExModel::restoreLastDeletedEvent(void)
{


    // Use constructor for speciation and extinction

    SpExBranchEvent* newEvent =
        new SpExBranchEvent(0.0, 0.0, 0.0, 0.0,
            _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _rng,
            _lastDeletedEventMapTime);

    newEvent->setLamInit(_lastDeletedEventLambdaInit);
    newEvent->setLamShift(_lastDeletedEventLambdaShift);
    newEvent->setMuInit(_lastDeletedEventMuInit);
    newEvent->setMuShift(_lastDeletedEventMuShift);



    // add the event to the branch history.
    //  ALWAYS done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->addEventToBranchHistory(newEvent);

    _eventCollection.insert(newEvent);

    // Event is now inserted into branch history:
    //  however, branch histories must be updated.

    forwardSetBranchHistories(newEvent);

    _tree->setMeanBranchSpeciation();
    _tree->setMeanBranchExtinction();

}



void SpExModel::changeNumberOfEventsMH(void)
{


    // Get old prior density of the data:
    double oldLogPrior = computeLogPrior();
    double newLogPrior = 0.0;
    int currState = (int)_eventCollection.size();
    int proposedState = 0;
    bool acceptMove = false;

    // Propose gains & losses equally if not on boundary (n = 0) events:

    // Current number of events on the tree, not counting root state:
    double K = (double)(_eventCollection.size());

    bool gain = _rng->uniformRv() <= 0.5;
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

        _tree->setMeanBranchSpeciation();
        _tree->setMeanBranchExtinction();

        double likBranches = computeLikelihoodBranches();
        double PropLnLik = likBranches;



#endif
        // Prior density on all parameters
        newLogPrior = computeLogPrior();


        // Prior ratio is eventRate / (k + 1)
        // but now, _eventCollection.size() == k + 1
        //  because event has already been added.
        // Here HR is just the prior ratio

        double logHR = log(_eventRate) - log(K + 1.0);

        // Now add log qratio

        logHR += log(qratio);

        double likeRatio = PropLnLik - getCurrLnLBranches();

        logHR += likeRatio;

        // Now for prior:
        logHR += (newLogPrior - oldLogPrior);

        // Now for jumping density of the bijection between parameter spaces:
        logHR -= _logQRatioJump;

        if (std::isinf(likeRatio) ) {

        } else
            acceptMove = acceptMetropolisHastings(logHR);




        if (acceptMove) {
            //std::cout << "gaining event in changeNumberOfEventsMH " << std::endl;
            //std::cout << "gainaccept" << computeLikelihoodBranches()  << std::endl;

            //addEventToTree();
            //std::cout << "Calliing isValid from ChangeNumberEvents::gain" << std::endl;

            bool isValidConfig = isEventConfigurationValid(_lastEventModified);

            if (isValidConfig) {
                // Update accept/reject statistics
                _acceptCount++;
                _acceptLast = 1;
            } else {
                // Need to get rid of event that was just gained...
                //std::cout << "Invalid event config from addEventToTree - deleting." << std::endl;
                deleteEventFromTree(_lastEventModified);
                _rejectCount++;
                _acceptLast = 0;

            }


        } else {
            //std::cout << "gainreject" << computeLikelihoodBranches() << std::endl;

            // Delete event.
            deleteEventFromTree(_lastEventModified);

            _rejectCount++;
            _acceptLast = 0;
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

        _tree->setMeanBranchSpeciation();
        _tree->setMeanBranchExtinction();

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
        double logHR = log(K) - log(_eventRate);

        // Now correct for proposal ratio:
        logHR += log(qratio);

        double likeRatio = PropLnLik - getCurrLnLBranches();


        logHR += likeRatio;

        // Now for prior:
        logHR += (newLogPrior - oldLogPrior);

        // Now for jumping density of the bijection between parameter spaces:
        logHR += _logQRatioJump;

        if (std::isinf(likeRatio) ) {

        } else
            acceptMove = acceptMetropolisHastings(logHR);


        //std::cout << "loss: " << acceptMove << "\t" << PropLnLik << "\tLT " << getCurrLnLTraits() + getCurrLnLBranches() << std::endl;

        if (acceptMove) {
            //std::cout << "loss accept, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;
            setCurrLnLTraits(likTraits);

            setCurrLnLBranches(likBranches);

            _acceptCount++;
            _acceptLast = 1;


        } else {
            //std::cout << "loss reject, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;

            restoreLastDeletedEvent();

            // speciation-extinction rates on branches automatically updated after restoreLastDeletedEvent()

            //std::cout << "loss reject restored, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << std::endl;
            _rejectCount++;
            _acceptLast = 0;
        }

    }

    //std::cout << currState << "\t" << proposedState << "\tAcc: " << acceptMove << "\tOP: " << oldLogPrior;
    //std::cout << "\tNP: " << newLogPrior << "\tqratio: " << _logQRatioJump << std::endl;

    incrementGeneration();

}

void SpExModel::moveEventMH(void)
{


    if (_eventCollection.size() > 0) {

        double localMoveProb = _localGlobalMoveRatio / (1 + _localGlobalMoveRatio);

        bool isLocalMove = _rng->uniformRv() <= localMoveProb;
        //std::cout << "is local: " << isLocalMove << std::endl;

        if (isLocalMove) {
            // Local move, with event drawn at random
            eventLocalMove();

#ifdef NO_DATA
            double likBranches = 0;
            double PropLnLik = likBranches;
#else
            _tree->setMeanBranchSpeciation();
            _tree->setMeanBranchExtinction();

            double likBranches = computeLikelihoodBranches();
            double PropLnLik = likBranches;



#endif

            double likeRatio = PropLnLik - getCurrLnLBranches();
            double logHR = likeRatio;

            // No longer protecting this as const bool

            bool acceptMove = false;
            bool isValid = false;
            //std::cout << "calling isValid from moveEventMH::local" << std::endl;
            isValid = isEventConfigurationValid(_lastEventModified);

            if (std::isinf(likeRatio) ) {

            } else if (isValid)
                acceptMove = acceptMetropolisHastings(logHR);
            else {
                //std::cout << "Invalid event configuration from LocalMove" << std::endl;
            }

            //const bool acceptMove = acceptMetropolisHastings(logHR);

            if (acceptMove == true) {

                setCurrLnLBranches(likBranches);
                _acceptCount++;
                _acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                _tree->setMeanBranchSpeciation();
                _tree->setMeanBranchExtinction();

                _rejectCount++;
                _acceptLast = 0;
            }



        } else {
            // global move, event drawn at random
            eventGlobalMove();
            //std::cout << "successful global move" << std::endl;
#ifdef NO_DATA
            double likBranches = 0;
            double PropLnLik = likBranches;
#else

            _tree->setMeanBranchSpeciation();
            _tree->setMeanBranchExtinction();

            double likBranches = computeLikelihoodBranches();
            double PropLnLik = likBranches;

#endif

            double likeRatio = PropLnLik - getCurrLnLBranches();
            double logHR = likeRatio;


            bool acceptMove = false;
            bool isValid = false;

            isValid = isEventConfigurationValid(_lastEventModified);

            if (std::isinf(likeRatio) ) {

            } else if (isValid)
                acceptMove = acceptMetropolisHastings(logHR);
            else {

            }

            if (acceptMove == true) {

                setCurrLnLBranches(likBranches);
                _acceptCount++;
                _acceptLast = 1;

            } else {
                // revert to previous state
                revertMovedEventToPrevious();

                _tree->setMeanBranchSpeciation();
                _tree->setMeanBranchExtinction();

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




void SpExModel::updateLambdaInitMH(void)
{

    //int n_events = _eventCollection.size() + 1;
    int toUpdate =_rng->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = static_cast<SpExBranchEvent*>(*myIt);
    }

    double oldRate = be->getLamInit();
    double cterm = exp( _updateLambdaInitScale * (_rng->uniformRv() - 0.5) );

    be->setLamInit(cterm * oldRate);

#ifdef NO_DATA
    double PropLnLik = 0;
#else

    _tree->setMeanBranchSpeciation();
    _tree->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

        
    double logPriorRatio = _prior->lambdaInitPrior(be->getLamInit());
    logPriorRatio -= _prior->lambdaInitPrior(oldRate);
    
    
    double LogProposalRatio = log(cterm);

    double likeRatio = PropLnLik - getCurrLnLBranches();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrLnLBranches(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
    } else {

        // revert to previous state
        be->setLamInit(oldRate);

        _tree->setMeanBranchSpeciation();
        _tree->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
    }

    incrementGeneration();


}

void SpExModel::updateLambdaShiftMH(void)
{

    //int n_events = _eventCollection.size() + 1;
    int toUpdate =_rng->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = static_cast<SpExBranchEvent*>(*myIt);
    }

    double oldLambdaShift = be->getLamShift();
    double newLambdaShift = oldLambdaShift +_rng->normalRv((double)0.0,
                            _updateLambdaShiftScale);
    be->setLamShift(newLambdaShift);


#ifdef NO_DATA
    double PropLnLik = 0;
#else

    _tree->setMeanBranchSpeciation();
    _tree->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

    double  logPriorRatio = _prior->lambdaShiftPrior(newLambdaShift);
    logPriorRatio -= _prior->lambdaShiftPrior(oldLambdaShift);
    
/*
    double  logPriorRatio =_rng->lnNormalPdf((double)0.0,
                            _settings->getLambdaShiftPrior(), newLambdaShift);
    logPriorRatio -=_rng->lnNormalPdf((double)0.0, _settings->getLambdaShiftPrior(),
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
        _acceptCount++;
        _acceptLast = 1;
    } else {

        // revert to previous state
        be->setLamShift(oldLambdaShift);

        _tree->setMeanBranchSpeciation();
        _tree->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
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
void SpExModel::updateTimeVariablePartitionsMH(void)
{

    //int n_events = _eventCollection.size() + 1;
    int toUpdate =_rng->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = static_cast<SpExBranchEvent*>(*myIt);
    } else {
        // event remains as root event-
    }

    if (be->getIsEventTimeVariable()) {




    } else if (!be->getIsEventTimeVariable()) {



    } else {
        // Should not be able to get here:
        std::cout << "Invalid _isEventTimeVariable in SpExModel::UpdateTimeVariablePartitionsMH"
             << std::endl;
        throw;
    }


}


void SpExModel::updateMuInitMH(void)
{

    //int n_events = _eventCollection.size() + 1;
    int toUpdate =_rng->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = static_cast<SpExBranchEvent*>(*myIt);
    }

    double oldRate = be->getMuInit();
    double cterm = exp( _updateMuInitScale * (_rng->uniformRv() - 0.5) );

    be->setMuInit(cterm * oldRate);

#ifdef NO_DATA
    double PropLnLik = 0;
#else

    _tree->setMeanBranchSpeciation();
    _tree->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

    double logPriorRatio = _prior->muInitPrior(be->getMuInit());
    logPriorRatio -= _prior->muInitPrior(oldRate);

    double LogProposalRatio = log(cterm);

    double likeRatio = PropLnLik - getCurrLnLBranches();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;
    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrLnLBranches(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
    } else {

        // revert to previous state
        be->setMuInit(oldRate);

        _tree->setMeanBranchSpeciation();
        _tree->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
    }

    incrementGeneration();


}


void SpExModel::updateMuShiftMH(void)
{

    //int n_events = _eventCollection.size() + 1;
    int toUpdate =_rng->sampleInteger(0, (int)_eventCollection.size());
    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);

    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
        for (int i = 1; i < toUpdate; i++)
            myIt++;

        be = static_cast<SpExBranchEvent*>(*myIt);
    }

    double oldMuShift = be->getMuShift();
    double newMuShift = oldMuShift +_rng->normalRv((double)0.0,
                        _updateMuShiftScale);

    be->setMuShift(newMuShift);

#ifdef NO_DATA
    double PropLnLik = 0;
#else

    _tree->setMeanBranchSpeciation();
    _tree->setMeanBranchExtinction();

    double PropLnLik = computeLikelihoodBranches();

#endif

    
    double logPriorRatio = _prior->muShiftPrior(newMuShift);
    logPriorRatio -= _prior->muShiftPrior(oldMuShift);
    
    

    double LogProposalRatio = 0.0;

    double likeRatio = PropLnLik - getCurrLnLBranches();

    double logHR = likeRatio +  logPriorRatio + LogProposalRatio;

    bool acceptMove = false;



    if (std::isinf(likeRatio) ) {

    } else
        acceptMove = acceptMetropolisHastings(logHR);


    if (acceptMove == true) {

        setCurrLnLBranches(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {

        // revert to previous state
        be->setMuShift(oldMuShift);

        _tree->setMeanBranchSpeciation();
        _tree->setMeanBranchExtinction();

        _acceptLast = 0;
        _rejectCount++;
    }

    incrementGeneration();


}




/*

 Metropolis-Hastings step to update Poisson event rate.
 Note that changing this rate does not affect the likelihood,
 so the priors and qratio determine acceptance rate.

 */
void SpExModel::updateEventRateMH(void)
{


    double oldEventRate = getEventRate();

    double cterm = exp( _updateEventRateScale * (_rng->uniformRv() - 0.5) );
    setEventRate(cterm * oldEventRate);

    
    double logPriorRatio = _prior->poissonRatePrior(getEventRate());
    logPriorRatio -= _prior->poissonRatePrior(oldEventRate);
    

    double logProposalRatio = log(cterm);


    // Experimental code:
    // Sample new event rate from prior directly with each step.
    //double newEventRate =_rng->exponentialRv(_poissonRatePrior);
    //setEventRate(newEventRate);
    //double LogPriorRatio = 0.0;
    //double logProposalRatio = 1.0;

    double logHR = logPriorRatio + logProposalRatio;
    const bool acceptMove = acceptMetropolisHastings(logHR);


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

}


double SpExModel::computeLikelihoodBranches(void)
{

    return computeLikelihoodBranchesByInterval();
}


double SpExModel::computeLikelihoodBranchesByInterval(void)
{


    double LnL = 0.0;


    if (_settings->getSampleFromPriorOnly())
        return 0.0;

    int numNodes = _tree->getNumberOfNodes();

    // LEft and right extinction probabilities for root node
    double rootEleft = 0.0;
    double rootEright = 0.0;


    for (int i = 0; i < numNodes; i++) {
        Node* xnode = _tree->getNodeFromDownpassSeq(i);
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
            if (xnode == _tree->getRoot())
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
            if _rng->uniformRv() <= 0.5){
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


            if (xnode == _tree->getRoot())
                rootEright = rEinit;
            // rDinit at this point should be FINAL value:
            // save as new variable, to keep clear:

            /* SHould be abele to ignore all of these calculations for the root node:
             Must also compute speciation rate for focal node. THis is a critical step.

             Since we are using approximations for the calculations on branches, we should set node speciation
             rate to be equivalent. Currently, I am not doing this - just computing exact rates
             at nodes.
            */

            if (xnode != _tree->getRoot()) {

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

double SpExModel::computeLogPrior(void)
{

    double logPrior = 0.0;

    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    logPrior += _prior->lambdaInitPrior(rootEvent->getLamInit());
    logPrior += _prior->lambdaShiftPrior(rootEvent->getLamShift());
    logPrior += _prior->muInitPrior(rootEvent->getMuInit());
    logPrior += _prior->muShiftPrior(rootEvent->getMuShift());
    
    int ctr = 0;

    for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
            i != _eventCollection.end(); i++) {

        SpExBranchEvent* event = static_cast<SpExBranchEvent*>(*i);

        logPrior += _prior->lambdaInitPrior(event->getLamInit());
        logPrior += _prior->lambdaShiftPrior(event->getLamShift());
        logPrior += _prior->muInitPrior(event->getMuInit());
        logPrior += _prior->muShiftPrior(event->getMuShift());
        
        ctr++;

    }

    // Here's prior density on the event rate:
    
    logPrior += _prior->poissonRatePrior(getEventRate());
    
    return logPrior;

}





bool SpExModel::acceptMetropolisHastings(const double lnR)
{
    const double r = safeExponentiation(SpExModel::mhColdness * lnR);
    return _rng->uniformRv() < r;
}




void SpExModel::initializeBranchHistories(Node* x)
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



void SpExModel::printStartAndEndEventStatesForBranch(Node* x)
{

    if (x != _tree->getRoot())
        std::cout << "Node: " << x << "\tAnc: " <<
             x->getBranchHistory()->getAncestralNodeEvent() << "\tevent: " <<
             x->getBranchHistory()->getNodeEvent() << std::endl;

    if (x->getLfDesc() != NULL)
        printStartAndEndEventStatesForBranch(x->getLfDesc());
    if (x->getRtDesc() != NULL)
        printStartAndEndEventStatesForBranch(x->getRtDesc());
}


void SpExModel::printBranchHistories(Node* x)
{

    if (x != _tree->getRoot()) {
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



double  SpExModel::getMHacceptanceRate(void)
{

    double arate = (double)_acceptCount / ((double)_acceptCount +
                                          (double)_rejectCount);

    return arate;

}


void  SpExModel::resetMHacceptanceParameters(void)
{
    _acceptCount = 0;
    _rejectCount = 0;

}



BranchEvent* SpExModel::getEventByIndex(int x)
{

    //int ctr = 0;
    std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
    for (int i = 0; i <= x; i++)
        myIt++;

    return (*myIt);
}


// adding log contrasts here.


void SpExModel::printExtinctionParams(void)
{

    if (_eventCollection.size() > 0) {
        for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
                i != _eventCollection.end(); i++) {
            SpExBranchEvent* event = static_cast<SpExBranchEvent*>(*i);
            std::cout << event << "\t" << event->getMuInit() << "\t" << event->getMuShift() << std::endl;
        }
    }

    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    std::cout << rootEvent << "\t" << rootEvent->getMuInit() << "\t" <<
        rootEvent->getMuShift() << std::endl << std::endl;
}


/*
 SpExModel::countTimeVaryingRatePartitions

    -counts number of time-varying rate partitions

 */
int SpExModel::countTimeVaryingRatePartitions(void)
{

    int count = 0;
    count += (int)_rootEvent->getIsEventTimeVariable();
    for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
            i != _eventCollection.end(); i++)
        count += (int)(*i)->getIsEventTimeVariable();
    return count;
}


/*
    Write event data to file for all events "on" tree
        at a given point in the MCMC chain


 */

void SpExModel::getEventDataString(std::stringstream& ss)
{

    ss << getGeneration() << ",";


    SpExBranchEvent* be = static_cast<SpExBranchEvent*>(_rootEvent);
    Node* xl = _tree->getRoot()->getRandomLeftTipNode();
    Node* xr = _tree->getRoot()->getRandomRightTipNode();
    ss << xl->getName() << "," << xr->getName() << "," << be->getAbsoluteTime() <<
       ",";
    ss << be->getLamInit() << "," << be->getLamShift() << "," << be->getMuInit() <<
       ",";
    ss << be->getMuShift();




    if (_eventCollection.size() > 0) {
        for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
                i != _eventCollection.end(); i++) {

            ss << "\n" << getGeneration() << ",";
            be = static_cast<SpExBranchEvent*>(*i);
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


bool SpExModel::isEventConfigurationValid(BranchEvent* be)
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
            std::cout << "problem in SpExModel::isEventConfigurationValid" << std::endl;
            exit(1);
        }


    }


    //std::cout << "leaving isEventConfigValid. Value: " << isValidConfig << std::endl;
    return isValidConfig;
}




double SpExModel::safeExponentiation(double x)
{
    if (x > 0.0)
        return 1.0;
    else if (x < -100.0)
        return 0.0;
    else
        return exp(x);
}


void SpExModel::debugLHcalculation(void)
{
    std::cout << "This does not currently support anything" << std::endl;

}









