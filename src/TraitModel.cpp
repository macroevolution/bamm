// defining this macro constrains analysis to NEGATIVE values
// for the beta shift parameter

//#define NEGATIVE_SHIFT_PARAM
#undef NEGATIVE_SHIFT_PARAM

#undef DEBUG  // This is a problem.


#include "TraitModel.h"
#include "Random.h"
#include "Settings.h"
#include "Tree.h"
#include "Node.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "TraitBranchEvent.h"
#include "BetaInitProposal.h"
#include "BetaShiftProposal.h"
#include "BetaTimeModeProposal.h"
#include "NodeStateProposal.h"
#include "Log.h"
#include "Prior.h"
#include "Stat.h"
#include "Tools.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <cstdlib>
#include <sstream>
#include <cmath>


TraitModel::TraitModel(Random& random, Settings& settings) :
    Model(random, settings)
{
#ifdef NEGATIVE_SHIFT_PARAM
    // Constrain beta shift to be zero or less than zero.
    if (_settings.getBetaShiftInit() > 0) {
        log(Error) << "Initial value of beta shift (betaShiftInit) cannot be\n"
            << "positive. This parameter is constrained to negative values\n";
        std::exit(1);
    }
#endif
    
    _sampleFromPriorOnly = _settings.get<bool>("sampleFromPriorOnly");

    double betaInit = _settings.get<double>("betaInit");
    double betaShiftInit = _settings.get<double>("betaShiftInit");

    bool isTimeVariable = false;
    double timeVarPrior = _settings.get<double>("betaIsTimeVariablePrior");
    if (timeVarPrior == 0.0) {
        isTimeVariable = false;
        if (betaShiftInit != 0.0) {
            exitWithError("betaShiftInit needs to be 0.0 if "
                "betaIsTimeVariablePrior is also 0.0");
        }
    } else if (timeVarPrior == 1.0) {
        isTimeVariable = true;
    } else {
        isTimeVariable = betaShiftInit != 0.0;
    }

    BranchEvent* x = new TraitBranchEvent(betaInit, betaShiftInit,
        isTimeVariable, _tree->getRoot(), _tree, _random, 0);
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the _rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event:
    forwardSetBranchHistories(_rootEvent);

    _tree->setMeanBranchTraitRates();

    // Initialize by previous event histories (or from initial event number)
    if (_settings.get<bool>("loadEventData")) {
        initializeModelFromEventDataFile(_settings.get("eventDataInfile"));
    } else {
        int initialNumberOfEvents = _settings.get<int>("initialNumberEvents");
        for (int i = 0; i < initialNumberOfEvents; i++) {
            addRandomEventToTree();
        }
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

    // Add proposals
    _proposals.push_back(new BetaInitProposal(random, settings, *this, _prior));
    _proposals.push_back
        (new BetaShiftProposal(random, settings, *this, _prior));
    _proposals.push_back(new NodeStateProposal(random, settings, *this));
    _proposals.push_back(new BetaTimeModeProposal(random, settings, *this));

 
    Model::calculateUpdateWeights();
 
    
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

    // TODO: Return true for now for time-variable
    return new TraitBranchEvent(betaInit, betaShift, true,
            x, _tree, _random, time);
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
    bool newIsTimeVariable = _prior.generateBetaIsTimeVariableFromPrior();

#ifdef NEGATIVE_SHIFT_PARAM
    newBetaShift = -fabs(newBetaShift);
    double dens_term = std::log(2.0);
#else
    double dens_term = 0.0;
#endif

    _logQRatioJump = 0.0;

    _logQRatioJump += _prior.betaInitPrior(newbeta);
    if (newIsTimeVariable) {
        _logQRatioJump += dens_term + _prior.betaShiftPrior(newBetaShift);
    }

    return new TraitBranchEvent(newbeta, newBetaShift, newIsTimeVariable,
        _tree->mapEventToTree(x), _tree, _random, x);
}


// TODO: test this : has not been checked
//     also not implemented in constructor for TraitModel yet
BranchEvent* TraitModel::newBranchEventWithParametersFromSettings(double x)
{
    
    // x is map time
    double newbeta = _settings.get<double>("betaInit0");
    double newBetaShift = _settings.get<double>("betaShift0");
    bool newIsTimeVariable = _prior.generateBetaIsTimeVariableFromPrior();

    
    // TODO: This needs to be refactored somewhere else
    // Computes the jump density for the addition of new parameters.
#ifdef NEGATIVE_SHIFT_PARAM
    newBetaShift = -fabs(newBetaShift);
    double dens_term = std::log(2.0);
#else
    double dens_term = 0.0;
#endif
    
    _logQRatioJump = 0.0;
    
    _logQRatioJump += _prior.betaInitPrior(newbeta);
    if (newIsTimeVariable) {
        _logQRatioJump += dens_term + _prior.betaShiftPrior(newBetaShift);
    }
    
    return new TraitBranchEvent(newbeta, newBetaShift, newIsTimeVariable,
                                _tree->mapEventToTree(x), _tree, _random, x);
    
}



void TraitModel::setDeletedEventParameters(BranchEvent* be)
{
    TraitBranchEvent* event = static_cast<TraitBranchEvent*>(be);

    _lastDeletedEventBetaInit = event->getBetaInit();
    _lastDeletedEventBetaShift = event->getBetaShift();
    _lastDeletedEventTimeVariable = event->isTimeVariable();
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
    return new TraitBranchEvent(_lastDeletedEventBetaInit,
        _lastDeletedEventBetaShift, _lastDeletedEventTimeVariable,
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _random,
        _lastDeletedEventMapTime);
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

    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        Node* xnode = postOrderNodes[i];
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
    double dens_term = std::log(2.0);
#else
    double dens_term = 0.0;
#endif

    double logPrior = 0.0;

    TraitBranchEvent* re = static_cast<TraitBranchEvent*>(_rootEvent);

    logPrior += _prior.betaInitRootPrior(re->getBetaInit());
    if (re->isTimeVariable()) {
        logPrior += dens_term + _prior.betaShiftRootPrior(re->getBetaShift());
    }

    for (std::set<BranchEvent*>::iterator i = _eventCollection.begin();
         i != _eventCollection.end(); ++i) {

        TraitBranchEvent* event = static_cast<TraitBranchEvent*>(*i);

        logPrior += _prior.betaInitPrior(event->getBetaInit());
        if (event->isTimeVariable()) {
            logPrior += dens_term +
                _prior.betaShiftPrior(event->getBetaShift());
        }

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
