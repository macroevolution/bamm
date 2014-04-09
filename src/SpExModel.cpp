#include <iostream>
#include <iomanip>
#include <cstdlib>
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


#define JUMP_VARIANCE_NORMAL 0.05


SpExModel::SpExModel(MbRandom* ranptr, Settings* sp) :
    Model(ranptr, sp),
    _lambdaInitProposal(*ranptr, *sp, *this, _prior),
    _lambdaShiftProposal(*ranptr, *sp, *this, _prior),
    _muInitProposal(*ranptr, *sp, *this, _prior),
    _muShiftProposal(*ranptr, *sp, *this, _prior)
{
    // Initialize tree
    if (_settings->get<bool>("useGlobalSamplingProbability")) {
        _tree->initializeSpeciationExtinctionModel
            (_settings->get<double>("globalSamplingFraction"));
    } else {
        _tree->initializeSpeciationExtinctionModel
            (_settings->get("sampleProbsFilename"));
    }

    _tree->setCanNodeHoldEventByDescCount
        (_settings->get<int>("minCladeSizeForShift"));
    _tree->setTreeMap(_tree->getRoot());

    // Initial values
    _lambdaInit0 = _settings->get<double>("lambdaInit0");
    _lambdaShift0 = _settings->get<double>("lambdaShift0");
    _muInit0 = _settings->get<double>("muInit0");
    _muShift0 = _settings->get<double>("muShift0");

    _sampleFromPriorOnly = _settings->get<bool>("sampleFromPriorOnly");

    // Parameter for splitting branch into pieces for numerical computation
    _segLength =
        _settings->get<double>("segLength") * _tree->maxRootToTipLength();
   
    BranchEvent* x =  new SpExBranchEvent
        (_lambdaInit0, _lambdaShift0, _muInit0, _muShift0,
            _tree->getRoot(), _tree, _rng, 0);
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the_rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event
    forwardSetBranchHistories(_rootEvent);

    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();

    // Initialize by previous event histories
    if (_settings->get<bool>("loadEventData")) {
        log() << "\nLoading model data from file.\n";
        initializeModelFromEventDataFile();
    }

    _extinctionProbMax = _settings->get<double>("extinctionProbMax");

    setCurrentLogLikelihood(computeLogLikelihood());

    if (std::isinf(getCurrentLogLikelihood())) {
        log(Error) << "Initial log-likelihood is infinity.\n"
            << "Please check your initial parameter values.\n";
        std::exit(1);
    }

    log() << "\nInitial log-likelihood: " << getCurrentLogLikelihood() << "\n";
    if (_sampleFromPriorOnly)
        log() << "Note that you have chosen to sample from prior only.\n";

    Model::finishConstruction();
}


void SpExModel::initializeSpecificUpdateWeights()
{
    _updateWeights.push_back(_settings->get<double>("updateRateLambda0"));
    _updateWeights.push_back(_settings->get<double>("updateRateLambdaShift"));
    _updateWeights.push_back(_settings->get<double>("updateRateMu0"));
    _updateWeights.push_back(_settings->get<double>("updateRateMuShift"));
}


Proposal* SpExModel::getSpecificProposal(int parameter)
{
    if (parameter == 3) {
        return &_lambdaInitProposal;
    } else if (parameter == 4) {
        return &_lambdaShiftProposal;
    } else if (parameter == 5) {
        return &_muInitProposal;
    } else if (parameter == 6) {
        return &_muShiftProposal;
    } else {
        // Should never get here
        log(Warning) << "Bad parameter to update.\n";
        return NULL;
    }
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
    double newLam = _prior.generateLambdaInitFromPrior();
    double newLambdaShift = _prior.generateLambdaShiftFromPrior();
    double newMu = _prior.generateMuInitFromPrior();
    double newMuShift = _prior.generateMuShiftFromPrior();
 
    // TODO: This needs to be refactored somewhere else
    // Computes the jump density for the addition of new parameters.
    _logQRatioJump = 0.0;    // Set to zero to clear previous values
    _logQRatioJump = _prior.lambdaInitPrior(newLam);
    _logQRatioJump += _prior.lambdaShiftPrior(newLambdaShift);
    _logQRatioJump += _prior.muInitPrior(newMu);
    _logQRatioJump += _prior.muShiftPrior(newMuShift);
    
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
    
    _logQRatioJump = _prior.lambdaInitPrior(_lastDeletedEventLambdaInit);
    _logQRatioJump += _prior.lambdaShiftPrior(_lastDeletedEventLambdaShift);
    _logQRatioJump += _prior.muInitPrior(_lastDeletedEventMuInit);
    _logQRatioJump += _prior.muShiftPrior(_lastDeletedEventMuShift);

    return _logQRatioJump;
}


BranchEvent* SpExModel::newBranchEventFromLastDeletedEvent()
{
    SpExBranchEvent* newEvent = new SpExBranchEvent(0.0, 0.0, 0.0, 0.0,
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _rng,
            _lastDeletedEventMapTime);

    newEvent->setLamInit(_lastDeletedEventLambdaInit);
    newEvent->setLamShift(_lastDeletedEventLambdaShift);
    newEvent->setMuInit(_lastDeletedEventMuInit);
    newEvent->setMuShift(_lastDeletedEventMuShift);

    return newEvent;
}


double SpExModel::computeLogLikelihood()
{
    return computeLogLikelihoodByInterval();
}



double SpExModel::computeLogLikelihoodByInterval()
{


    double LnL = 0.0;


    if (_sampleFromPriorOnly)
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

            if (lEinit > _extinctionProbMax || rEinit > _extinctionProbMax) {
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




double SpExModel::computeLogPrior()
{
    double logPrior = 0.0;

    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    logPrior += _prior.lambdaInitRootPrior(rootEvent->getLamInit());
    logPrior += _prior.lambdaShiftRootPrior(rootEvent->getLamShift());
    logPrior += _prior.muInitRootPrior(rootEvent->getMuInit());
    logPrior += _prior.muShiftRootPrior(rootEvent->getMuShift());

    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        SpExBranchEvent* event = static_cast<SpExBranchEvent*>(*it);

        logPrior += _prior.lambdaInitPrior(event->getLamInit());
        logPrior += _prior.lambdaShiftPrior(event->getLamShift());
        logPrior += _prior.muInitPrior(event->getMuInit());
        logPrior += _prior.muShiftPrior(event->getMuShift());
    }

    // Here's prior density on the event rate
    logPrior += _prior.poissonRatePrior(_eventRate);
    
    return logPrior;
}
