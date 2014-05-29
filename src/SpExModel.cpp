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
#include "Random.h"
#include "Node.h"
#include "Tree.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "SpExBranchEvent.h"
#include "Settings.h"
#include "Log.h"
#include "Prior.h"
#include "Tools.h"

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


SpExModel::SpExModel(Random& random, Settings& settings) :
    Model(random, settings),
    _lambdaInitProposal(random, settings, *this, _prior),
    _lambdaShiftProposal(random, settings, *this, _prior),
    _muInitProposal(random, settings, *this, _prior),
    _muShiftProposal(random, settings, *this, _prior)
{
    // Initial values
    _lambdaInit0 = _settings.get<double>("lambdaInit0");
    _lambdaShift0 = _settings.get<double>("lambdaShift0");
    _muInit0 = _settings.get<double>("muInit0");
    _muShift0 = _settings.get<double>("muShift0");

    _sampleFromPriorOnly = _settings.get<bool>("sampleFromPriorOnly");

    // Parameter for splitting branch into pieces for numerical computation
    _segLength =
        _settings.get<double>("segLength") * _tree->maxRootToTipLength();
   
    BranchEvent* x =  new SpExBranchEvent
        (_lambdaInit0, _lambdaShift0, _muInit0, _muShift0,
            _tree->getRoot(), _tree, _random, 0);
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the_rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event
    forwardSetBranchHistories(_rootEvent);

    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();

    // Initialize by previous event histories
    if (_settings.get<bool>("loadEventData")) {
        initializeModelFromEventDataFile(_settings.get("eventDataInfile"));
    }

    _extinctionProbMax = _settings.get<double>("extinctionProbMax");

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
    _updateWeights.push_back(_settings.get<double>("updateRateLambda0"));
    _updateWeights.push_back(_settings.get<double>("updateRateLambdaShift"));
    _updateWeights.push_back(_settings.get<double>("updateRateMu0"));
    _updateWeights.push_back(_settings.get<double>("updateRateMuShift"));
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


void SpExModel::setRootEventWithReadParameters
    (const std::vector<std::string>& parameters)
{
    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    rootEvent->setLamInit(lambdaInitParameter(parameters));
    rootEvent->setLamShift(lambdaShiftParameter(parameters));
    rootEvent->setMuInit(muInitParameter(parameters));
    rootEvent->setMuShift(muShiftParameter(parameters));
}


BranchEvent* SpExModel::newBranchEventWithReadParameters
    (Node* x, double time, const std::vector<std::string>& parameters)
{
    double lambdaInit = lambdaInitParameter(parameters);
    double lambdaShift = lambdaShiftParameter(parameters);
    double muInit = muInitParameter(parameters);
    double muShift = muShiftParameter(parameters);

    return new SpExBranchEvent(lambdaInit, lambdaShift,
        muInit, muShift, x, _tree, _random, time);
}


double SpExModel::lambdaInitParameter
    (const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[0]);
}


double SpExModel::lambdaShiftParameter
    (const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[1]);
}


double SpExModel::muInitParameter(const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[2]);
}


double SpExModel::muShiftParameter(const std::vector<std::string>& parameters)
{
    return convert_string<double>(parameters[3]);
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
        newMuShift, _tree->mapEventToTree(x), _tree, _random, x);
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
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _random,
            _lastDeletedEventMapTime);

    newEvent->setLamInit(_lastDeletedEventLambdaInit);
    newEvent->setLamShift(_lastDeletedEventLambdaShift);
    newEvent->setMuInit(_lastDeletedEventMuInit);
    newEvent->setMuShift(_lastDeletedEventMuShift);

    return newEvent;
}


double SpExModel::computeLogLikelihood()
{
    if (_sampleFromPriorOnly)
        return 0.0;

    double logLikelihood = 0.0;

    int numNodes = _tree->getNumberOfNodes();

    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        Node* node = postOrderNodes[i];
        if (node->isInternal()) {
            logLikelihood += computeSpExProbBranch(node->getLfDesc());
            logLikelihood += computeSpExProbBranch(node->getRtDesc());

            // Does not include root node, so it is conditioned
            // on basal speciation event occurring:
            if (node != _tree->getRoot()) {
                logLikelihood  += log(node->getNodeLambda());
                node->setDinit(1.0);
            }
        }
    }

    return logLikelihood;
}


double SpExModel::computeSpExProbBranch(Node* node)
{
    double logLikelihood = 0.0;

    double D0 = node->getDinit();    // Initial speciation probability
    double E0 = node->getEinit();    // Initial extinction probability

    double startTime = node->getBrlen();
    double endTime = node->getBrlen();

    while (startTime > 0) {
        startTime -= _segLength;
        if (startTime < 0) {
            startTime = 0.0;
        }

        double deltaT = endTime - startTime;

        double curLam = node->computeSpeciationRateIntervalRelativeTime
            (startTime, endTime);

        double curMu = node->computeExtinctionRateIntervalRelativeTime
            (startTime, endTime);

        double spProb = 0.0;
        double exProb = 0.0;

        // Compute speciation and extinction probabilities and store them
        // in spProb and exProb (through reference passing)
        computeSpExProb(spProb, exProb, curLam, curMu, D0, E0, deltaT);

        logLikelihood += std::log(spProb);

        D0 = 1.0;
        E0 = exProb;

        endTime = startTime;
    }

    // What to use as E0 for the start of next downstream branch?
    // At this point, E0 is actually the E0 for the parent node.
    // Here, we will arbitrarily take this to be the left extinction value
    Node* parent = node->getAnc();
    if (node == parent->getLfDesc()) {
        parent->setEinit(E0);
    }

    // Clearly a problem if extinction values approaching/equaling 1.0
    // If so, set to -Inf, leading to automatic rejection of state
    if (E0 > _extinctionProbMax) {
        return -INFINITY;
    }

    // To CONDITION on survival, uncomment the line below;
    // or to NOT condition on survival, comment the line below
    if (parent == _tree->getRoot()) {
        logLikelihood -= std::log(1.0 - E0);
    }

    return logLikelihood;
}


// If we let
//
//     M = mu - lam
//     L = lam * (1.0 - E0)
//     E = e^(M * deltaT)
//     m = E * M
//     d = L * (1 - E) - m
//
// then the speciation equation from PLoS paper
//
//                    e^((mu - lam) * deltaT) * D0 * (lam - u)^2
//     D(t) = --------------------------------------------------------------
//            [lam - lam * E0 + e^((mu - lam) * deltaT) * (E0 * lam - mu)]^2
//
// becomes
//
//            D0 * m
//     D(t) = ------
//              d^2
//
// and the extinction equation
//
//                               (1 - E0) * (lam - mu)
//     E(t) = 1 - ----------------------------------------------------------
//                (1 - E0) * lam - e^(-(lam - mu) * deltaT) * (mu - lam * E0)
//
// becomes
//
//                (1 - E0) * M
//     E(t) = 1 + ------------
//                      d

void SpExModel::computeSpExProb(double& spProb, double& exProb,
    double lambda, double mu, double D0, double E0, double deltaT)
{
    double M = mu - lambda;
    double L = lambda * (1.0 - E0);
    double E = std::exp(M * deltaT);
    double m = E * M;
    double d = L * (1.0 - E) - m;

    spProb = (D0 * m * M) / sqr(d);
    exProb = 1.0 + (1.0 - E0) * M / d;
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
