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
    _lastEventModified = x;

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
        initializeModelFromEventDataFile();
    }

    setCurrentLogLikelihood(computeLogLikelihood());

    log() << "\nInitial log-likelihood: " << getCurrentLogLikelihood() << "\n";
    if (_settings->getSampleFromPriorOnly()) {
        log() << "Note that you have chosen to sample from prior only.\n";
    }
}


void TraitModel::readModelSpecificParameters(std::ifstream &inputFile)
{
    inputFile >> _readBetaInit;
    inputFile >> _readBetaShift;
}


void TraitModel::setRootEventWithReadParameters()
{
    TraitBranchEvent* rootEvent = static_cast<TraitBranchEvent*>(_rootEvent);

    rootEvent->setBetaInit(_readBetaInit);
    rootEvent->setBetaShift(_readBetaShift);
}


BranchEvent* TraitModel::newBranchEventWithReadParameters(Node* x, double time)
{
    return new TraitBranchEvent(_readBetaInit, _readBetaShift,
        x, _tree, _rng, time);
}


void TraitModel::setMeanBranchParameters()
{
    _tree->setMeanBranchTraitRates();
}


BranchEvent* TraitModel::newBranchEventWithRandomParameters(double x)
{
    // Sample beta and beta shift from prior
    double newbeta = _prior->generateBetaInitFromPrior();
    double newBetaShift = _prior->generateBetaShiftFromPrior();
    
#ifdef NEGATIVE_SHIFT_PARAM
    newBetaShift = -fabs(newBetaShift);
    double dens_term = log(2.0);
#else
    double dens_term = 0.0;
#endif
    
    _logQRatioJump = 0.0;
    
    _logQRatioJump += _prior->betaInitPrior(newbeta);
    _logQRatioJump += dens_term + _prior->betaShiftPrior(newBetaShift);
    
    return new TraitBranchEvent(newbeta, newBetaShift,
        _tree->mapEventToTree(x), _tree, _rng, x);
}


void TraitModel::setDeletedEventParameters(BranchEvent* be)
{
    TraitBranchEvent* event = static_cast<TraitBranchEvent*>(be);

    _lastDeletedEventBetaInit = event->getBetaInit();
    _lastDeletedEventBetaShift = event->getBetaShift();
}


double TraitModel::calculateLogQRatioJump()
{
    double _logQRatioJump = 0.0;
    
    _logQRatioJump = _prior->betaInitPrior(_lastDeletedEventBetaInit);
    _logQRatioJump += _prior->betaShiftPrior(_lastDeletedEventBetaShift);

    return _logQRatioJump;
}
    

BranchEvent* TraitModel::newBranchEventFromLastDeletedEvent()
{
    TraitBranchEvent* newEvent = new TraitBranchEvent(0.0, 0.0,
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _rng,
            _lastDeletedEventMapTime);

    newEvent->setBetaInit(_lastDeletedEventBetaInit);
    newEvent->setBetaShift(_lastDeletedEventBetaShift);

    return newEvent;
}


void TraitModel::updateBetaMH(void)
{
    int toUpdate = _rng->sampleInteger(0, (int)_eventCollection.size());
    
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);
    
    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
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
    
    double PropLnLik = computeLogLikelihood();
    
#endif

    double logPriorRatio = 0.0;
    if (be == _rootEvent){
        logPriorRatio = _prior->betaInitRootPrior(be->getBetaInit());
        logPriorRatio -= _prior->betaInitRootPrior(oldRate);
    } else {
        logPriorRatio = _prior->betaInitPrior(be->getBetaInit());
        logPriorRatio -= _prior->betaInitPrior(oldRate);
    }

    double LogProposalRatio = log(cterm);
    
    double likeRatio = PropLnLik - getCurrentLogLikelihood();
    
    double logHR = likeRatio + logPriorRatio + LogProposalRatio;    

    bool acceptMove = acceptMetropolisHastings(logHR);
    
    //  std::cout << getGeneration() << "\tL1: " << startLH << "\tL2: " << getCurrentLogLikelihood() << std::endl;
    
    
    if (acceptMove == true) {
        //std::cout << "accept: " << oldRate << "\t" << be->getBetaInit() << std::endl;
        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;
        
    } else {
        
        be->setBetaInit(oldRate);
        _tree->setMeanBranchTraitRates();
        _acceptLast = 0;
        _rejectCount++;
        
    }
    
    /*if (!acceptMove){
     std::cout << std::endl;
     std::cout << startLL << "\tCurr: " << getCurrentLogLikelihood() << "\tcalc: " << computeLogLikelihood() << std::endl;
     }*/
    
    incrementGeneration();
    
}


void TraitModel::updateBetaShiftMH(void)
{
    int toUpdate = _rng->sampleInteger(0, (int)_eventCollection.size());
    
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(_rootEvent);
    
    if (toUpdate > 0) {
        std::set<BranchEvent*>::iterator myIt = _eventCollection.begin();
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
    double PropLnLik = computeLogLikelihood();
    
#endif

    double logPriorRatio = 0.0;
    if (be == _rootEvent){
        logPriorRatio = _prior->betaShiftRootPrior(be->getBetaShift());
        logPriorRatio -= _prior->betaShiftRootPrior(oldShift);
    } else {
        logPriorRatio = _prior->betaShiftPrior(be->getBetaShift());
        logPriorRatio -= _prior->betaShiftPrior(oldShift);
    }

    double LogProposalRatio = 0.0;
    
    double likeRatio = PropLnLik - getCurrentLogLikelihood();
    
    double logHR = likeRatio + logPriorRatio + LogProposalRatio;
    
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    if (acceptMove == true) {
        
        setCurrentLogLikelihood(PropLnLik);
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
    double PropLnLik = getCurrentLogLikelihood() - oldTriadLogLik + newTriadLoglik;

#endif

    // set to zero for now...flat prior, unifor (see below).
    double LogPriorRatio = 0.0;


    double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0

    double likeRatio = PropLnLik - getCurrentLogLikelihood();


    double logHR = LogPriorRatio + logProposalRatio + likeRatio;
    bool acceptMove = acceptMetropolisHastings(logHR);

    // Here we do prior calculation to avoid computing infinity...
    if (newstate > _settings->getTraitPriorMax() ||
            newstate < _settings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrentLogLikelihood(PropLnLik);
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
    double PropLnLik = getCurrentLogLikelihood() - oldTriadLogLik + newTriadLoglik;



#endif

    // set to zero for now...flat prior, unifor (see below).
    double LogPriorRatio = 0.0;


    double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0

    double likeRatio = PropLnLik - getCurrentLogLikelihood();


    double logHR = LogPriorRatio + logProposalRatio + likeRatio;
    bool acceptMove = acceptMetropolisHastings(logHR);

    // Here we do prior calculation to avoid computing infinity...
    if (newstate > _settings->getTraitPriorMax() ||
            newstate < _settings->getTraitPriorMin())
        acceptMove = false;
    if (acceptMove == true) {
        //continue
        setCurrentLogLikelihood(PropLnLik);
        _acceptCount++;
        _acceptLast = 1;

    } else {
        xnode->setTraitValue(oldstate);

        _rejectCount++;
        _acceptLast = 0;
    }

    incrementGeneration();





}


double TraitModel::computeLogLikelihood()
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
    
    logPrior += _prior->betaInitRootPrior(re->getBetaInit());
    logPrior += dens_term + _prior->betaShiftRootPrior(re->getBetaShift());

    for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
         i != _eventCollection.end(); ++i) {

        TraitBranchEvent* event = static_cast<TraitBranchEvent*>(*i);
        
        logPrior += _prior->betaInitPrior(event->getBetaInit());
        logPrior += dens_term + _prior->betaShiftPrior(event->getBetaShift());
        
    }
    
    // and prior on number of events:
    
    logPrior += _prior->poissonRatePrior(getEventRate());
    
    return logPrior;
    
}


void TraitModel::getSpecificEventDataString
    (std::stringstream& ss, BranchEvent* event)
{
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(event);

    ss << be->getBetaInit() << ","
    << be->getBetaShift();
}


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
