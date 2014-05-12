// defining this macro constrains analysis to NEGATIVE values
// for the beta shift parameter

//#define NEGATIVE_SHIFT_PARAM
#undef NEGATIVE_SHIFT_PARAM

#undef DEBUG  // This is a problem.


#include "TraitModel.h"


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>
#include <algorithm>

#include "Random.h"
#include "Node.h"
#include "Tree.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "TraitBranchEvent.h"
#include "Settings.h"
#include "Log.h"
#include "Prior.h"
#include "Stat.h"
#include "Tools.h"


TraitModel::TraitModel(Random& random, Settings* settings) :
    Model(random, settings),
    _betaInitProposal(random, *settings, *this, _prior),
    _betaShiftProposal(random, *settings, *this, _prior),
    _nodeStateProposal(random, *settings, *this)
{
#ifdef NEGATIVE_SHIFT_PARAM
    // Constrain beta shift to be zero or less than zero.
    if (_settings->getBetaShiftInit() > 0) {
        log(Error) << "Initial value of beta shift (betaShiftInit) cannot be\n"
            << "positive. This parameter is constrained to negative values\n";
        std::exit(1);
    }
#endif
    
    _sampleFromPriorOnly = _settings->get<bool>("sampleFromPriorOnly");

    BranchEvent* x = new TraitBranchEvent(_settings->get<double>("betaInit"),
        _settings->get<double>("betaShiftInit"),
        _tree->getRoot(), _tree, _random, 0);
    _rootEvent = x;
    _lastEventModified = x;

    TraitBranchEvent* traitRootEvent =
        static_cast<TraitBranchEvent*>(_rootEvent);

    log() << "\nRoot beta: " << traitRootEvent->getBetaInit() << "\t"
          << _settings->get<double>("betaInit") << "\t"
          << "Shift: " << traitRootEvent->getBetaShift() << "\n";

    // Set NodeEvent of root node equal to the _rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event:
    forwardSetBranchHistories(_rootEvent);

    _tree->setMeanBranchTraitRates();

    if (_settings->get<bool>("loadEventData")) {
        initializeModelFromEventDataFile(_settings->get("eventDataInfile"));
    }

    setCurrentLogLikelihood(computeLogLikelihood());

    // TODO: Code duplication with SpExModel
    if (std::isinf(getCurrentLogLikelihood())) {
        log(Error) << "Initial log-likelihood is infinity.\n"
            << "Please check your initial parameter values.\n";
        std::exit(1);
    }

    log() << "\nInitial log-likelihood: " << getCurrentLogLikelihood() << "\n";
    if (_sampleFromPriorOnly) {
        log() << "Note that you have chosen to sample from prior only.\n";
    }

    Model::finishConstruction();
}


void TraitModel::initializeSpecificUpdateWeights()
{
    _updateWeights.push_back(_settings->get<double>("updateRateBeta0"));
    _updateWeights.push_back(_settings->get<double>("updateRateBetaShift"));
    _updateWeights.push_back(_settings->get<double>("updateRateNodeState"));
}


Proposal* TraitModel::getSpecificProposal(int parameter)
{
    if (parameter == 3) {
        return &_betaInitProposal;
    } else if (parameter == 4) {
        return &_betaShiftProposal;
    } else if (parameter == 5) {
        return &_nodeStateProposal;
    } else {
        // Should never get here
        log(Warning) << "Bad parameter to update.\n";
        return NULL;
    }
}


void TraitModel::setRootEventWithReadParameters
    (const std::vector<std::string>& parameters)
{
    TraitBranchEvent* rootEvent = static_cast<TraitBranchEvent*>(_rootEvent);

    rootEvent->setBetaInit(betaInitParameter(parameters));
    rootEvent->setBetaShift(betaShiftParameter(parameters));
}


BranchEvent* TraitModel::newBranchEventWithReadParameters
    (Node* x, double time, const std::vector<std::string>& parameters)
{
    double betaInit = betaInitParameter(parameters);
    double betaShift = betaShiftParameter(parameters);

    return new TraitBranchEvent(betaInit, betaShift, x, _tree, _random, time);
}


double TraitModel::betaInitParameter(const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[0]);
}


double TraitModel::betaShiftParameter
    (const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[1]);
}


void TraitModel::setMeanBranchParameters()
{
    _tree->setMeanBranchTraitRates();
}


BranchEvent* TraitModel::newBranchEventWithRandomParameters(double x)
{
    // Sample beta and beta shift from prior
    double newbeta = _prior.generateBetaInitFromPrior();
    double newBetaShift = _prior.generateBetaShiftFromPrior();
    
#ifdef NEGATIVE_SHIFT_PARAM
    newBetaShift = -fabs(newBetaShift);
    double dens_term = log(2.0);
#else
    double dens_term = 0.0;
#endif
    
    _logQRatioJump = 0.0;
    
    _logQRatioJump += _prior.betaInitPrior(newbeta);
    _logQRatioJump += dens_term + _prior.betaShiftPrior(newBetaShift);
    
    return new TraitBranchEvent(newbeta, newBetaShift,
        _tree->mapEventToTree(x), _tree, _random, x);
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
    
    _logQRatioJump = _prior.betaInitPrior(_lastDeletedEventBetaInit);
    _logQRatioJump += _prior.betaShiftPrior(_lastDeletedEventBetaShift);

    return _logQRatioJump;
}
    

BranchEvent* TraitModel::newBranchEventFromLastDeletedEvent()
{
    TraitBranchEvent* newEvent = new TraitBranchEvent(0.0, 0.0,
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _random,
            _lastDeletedEventMapTime);

    newEvent->setBetaInit(_lastDeletedEventBetaInit);
    newEvent->setBetaShift(_lastDeletedEventBetaShift);

    return newEvent;
}


double TraitModel::computeLogLikelihood()
{

    double LnL = 0.0;

    //Node * tmpnode = _tree->getRoot()->getLfDesc();

    if (_sampleFromPriorOnly)
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

            LnL += Stat::lnNormalPDF(delta, 0.0, std::sqrt(var));

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


    if (_sampleFromPriorOnly)
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
            double var = x->getLfDesc()->getBrlen() *
                x->getLfDesc()->getMeanBeta();
            logL += Stat::lnNormalPDF(delta, 0.0, std::sqrt(var));
        }


        if (x->getRtDesc()->getCanHoldEvent() == true) {
            // computation for right descendant branch
            double delta = x->getRtDesc()->getTraitValue() - x->getTraitValue();
            double var = x->getRtDesc()->getBrlen() *
                x->getRtDesc()->getMeanBeta();
            logL += Stat::lnNormalPDF(delta, 0.0, std::sqrt(var));
        }


        // computation for ancestral branch (unless == root)

        if (x != _tree->getRoot()) {

            double delta = x->getTraitValue() - x->getAnc()->getTraitValue();
            double var = x->getBrlen() * x->getMeanBeta();
            logL += Stat::lnNormalPDF(delta, 0.0, std::sqrt(var));
        }
    }

#ifdef DEBUG
    std::cout << "Leaving computeTriadLikelihood: Node : " << x << std::endl;
#endif

    return logL;

}


double TraitModel::computeLogPrior()
{
#ifdef NEGATIVE_SHIFT_PARAM
    double dens_term = log(2.0);
#else
    double dens_term = 0.0;
#endif
    
    double logPrior = 0.0;

    TraitBranchEvent* re = static_cast<TraitBranchEvent*>(_rootEvent);
    
    logPrior += _prior.betaInitRootPrior(re->getBetaInit());
    logPrior += dens_term + _prior.betaShiftRootPrior(re->getBetaShift());

    for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
         i != _eventCollection.end(); ++i) {

        TraitBranchEvent* event = static_cast<TraitBranchEvent*>(*i);
        
        logPrior += _prior.betaInitPrior(event->getBetaInit());
        logPrior += dens_term + _prior.betaShiftPrior(event->getBetaShift());
        
    }
    
    // and prior on number of events:
    
    logPrior += _prior.poissonRatePrior(getEventRate());
    
    return logPrior;
    
}


void TraitModel::getSpecificEventDataString
    (std::stringstream& ss, BranchEvent* event)
{
    TraitBranchEvent* be = static_cast<TraitBranchEvent*>(event);

    ss << be->getBetaInit() << ","
    << be->getBetaShift();
}
