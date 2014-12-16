#include "SpExModel.h"
#include "Model.h"
#include "Random.h"
#include "Settings.h"
#include "Tree.h"
#include "Node.h"
#include "BranchHistory.h"
#include "BranchEvent.h"
#include "SpExBranchEvent.h"
#include "LambdaInitProposal.h"
#include "LambdaShiftProposal.h"
#include "MuInitProposal.h"
#include "MuShiftProposal.h"
#include "LambdaTimeModeProposal.h"
#include "PreservationRateProposal.h"

#include "Log.h"
#include "Prior.h"
#include "Tools.h"

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
   

#define JUMP_VARIANCE_NORMAL 0.05


SpExModel::SpExModel(Random& random, Settings& settings) :
    Model(random, settings)
{
    // Initial values
    _lambdaInit0 = _settings.get<double>("lambdaInit0");
    _lambdaShift0 = _settings.get<double>("lambdaShift0");
    _muInit0 = _settings.get<double>("muInit0");
    _muShift0 = _settings.get<double>("muShift0");

    
    // Initialize fossil preservation rate:
    _preservationRate = _settings.get<double>("preservationRateInit");
    
    
    
    // Observation time of tree:
    _observationTime = _settings.get<double>("observationTime");
    if (_observationTime <= 0){
        _observationTime = _tree->getAge();
    }else if ( _observationTime < _tree->getAge() ){
        std::cout << "WARNING: invalid initial observation time" << std::endl;
        std::cout << "\t... setting observationTime to tree MAX TIME" << std::endl;
        _observationTime = _tree->getAge();
    }
    
    std::cout << "Observation time (from root): " << _observationTime << std::endl;
    
    // initialize occurrence count:
    _numberOccurrences = _settings.get<double>("numberOccurrences");
    if (_numberOccurrences > 0){
        // Good. User-defined.
    }else{
        // count fossil occurrences in tree.
        // Must be able to handle extant-only tree, eg zero value
        std::cout << "Must have valid number of occurrences specified" << std::endl;
        throw;
    }
    
    
    
    //_tree->printInitialSpeciationExtinctionRates();
     
    double timeVarPrior = _settings.get<double>("lambdaIsTimeVariablePrior");
    if (timeVarPrior == 0.0) {
        _initialLambdaIsTimeVariable = false;
        if (_lambdaShift0 != 0.0) {
            exitWithError("lambdaShift0 needs to be 0.0 if "
                "lambdaIsTimeVariablePrior is also 0.0");
        }
    } else if (timeVarPrior == 1.0) {
        _initialLambdaIsTimeVariable = true;
    } else {
        _initialLambdaIsTimeVariable = _lambdaShift0 != 0.0;
    }

    _sampleFromPriorOnly = _settings.get<bool>("sampleFromPriorOnly");

    // Parameter for splitting branch into pieces for numerical computation
    _segLength =
        _settings.get<double>("segLength") * _tree->maxRootToTipLength();

    BranchEvent* x =  new SpExBranchEvent(_lambdaInit0, _lambdaShift0,
        _muInit0, _muShift0, _initialLambdaIsTimeVariable,
        _tree->getRoot(), _tree, _random, 0);
    _rootEvent = x;
    _lastEventModified = x;

    // Set NodeEvent of root node equal to the_rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event
    forwardSetBranchHistories(_rootEvent);

    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();

    // Initialize by previous event histories (or from initial event number)
    if (_settings.get<bool>("loadEventData")) {
        initializeModelFromEventDataFile(_settings.get("eventDataInfile"));
    } else {
        int initialNumberOfEvents = _settings.get<int>("initialNumberEvents");
        for (int i = 0; i < initialNumberOfEvents; i++) {
            addRandomEventToTree();
        }
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

    // Add proposals
    _proposals.push_back
        (new LambdaInitProposal(random, settings, *this, _prior));
    _proposals.push_back
        (new LambdaShiftProposal(random, settings, *this, _prior));
    _proposals.push_back(new MuInitProposal(random, settings, *this, _prior));
    _proposals.push_back(new MuShiftProposal(random, settings, *this, _prior));
    _proposals.push_back(new LambdaTimeModeProposal(random, settings, *this));

    _proposals.push_back(new PreservationRateProposal(random, settings, *this, _prior));

    Model::calculateUpdateWeights();
}


void SpExModel::setRootEventWithReadParameters
    (const std::vector<std::string>& parameters)
{
    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    rootEvent->setLamInit(lambdaInitParameter(parameters));
    rootEvent->setLamShift(lambdaShiftParameter(parameters));
    rootEvent->setMuInit(muInitParameter(parameters));
    rootEvent->setMuShift(muShiftParameter(parameters));
    // TODO: Add time variable/constant
}


BranchEvent* SpExModel::newBranchEventWithReadParameters
    (Node* x, double time, const std::vector<std::string>& parameters)
{
    double lambdaInit = lambdaInitParameter(parameters);
    double lambdaShift = lambdaShiftParameter(parameters);
    double muInit = muInitParameter(parameters);
    double muShift = muShiftParameter(parameters);

    // TODO: Fix reading of parameters (for now, send true for time-variable)
    return new SpExBranchEvent(lambdaInit, lambdaShift,
        muInit, muShift, true, x, _tree, _random, time);
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
    double newMu = _prior.generateMuInitFromPrior();
    double newMuShift = _prior.generateMuShiftFromPrior();
    bool newIsTimeVariable = _prior.generateLambdaIsTimeVariableFromPrior();

    double newLambdaShift = 0.0;
    if (newIsTimeVariable) {
        newLambdaShift = _prior.generateLambdaShiftFromPrior();
    }

    // TODO: This needs to be refactored somewhere else
    // Computes the jump density for the addition of new parameters.
    _logQRatioJump = 0.0;    // Set to zero to clear previous values
    _logQRatioJump = _prior.lambdaInitPrior(newLam);
    if (newIsTimeVariable) {
        _logQRatioJump += _prior.lambdaShiftPrior(newLambdaShift);
    }
    _logQRatioJump += _prior.muInitPrior(newMu);
    _logQRatioJump += _prior.muShiftPrior(newMuShift);

    return new SpExBranchEvent(newLam, newLambdaShift, newMu,
        newMuShift, newIsTimeVariable, _tree->mapEventToTree(x),
        _tree, _random, x);
}


void SpExModel::setDeletedEventParameters(BranchEvent* be)
{
    SpExBranchEvent* event = static_cast<SpExBranchEvent*>(be);

    _lastDeletedEventLambdaInit = event->getLamInit();
    _lastDeletedEventLambdaShift = event->getLamShift();
    _lastDeletedEventMuInit = event->getMuInit();
    _lastDeletedEventMuShift = event->getMuShift();
    _lastDeletedEventTimeVariable = event->isTimeVariable();
}


double SpExModel::calculateLogQRatioJump()
{
    double _logQRatioJump = 0.0;

    _logQRatioJump = _prior.lambdaInitPrior(_lastDeletedEventLambdaInit);
    // TODO: Should something be checked for time variable/constant?
    _logQRatioJump += _prior.lambdaShiftPrior(_lastDeletedEventLambdaShift);
    _logQRatioJump += _prior.muInitPrior(_lastDeletedEventMuInit);
    _logQRatioJump += _prior.muShiftPrior(_lastDeletedEventMuShift);

    return _logQRatioJump;
}


BranchEvent* SpExModel::newBranchEventFromLastDeletedEvent()
{
    return new SpExBranchEvent(_lastDeletedEventLambdaInit,
        _lastDeletedEventLambdaShift, _lastDeletedEventMuInit,
        _lastDeletedEventMuShift, _lastDeletedEventTimeVariable,
        _tree->mapEventToTree(_lastDeletedEventMapTime), _tree, _random,
        _lastDeletedEventMapTime);
}


// TODO: Not transparent, but this is where
//  Di for internal nodes is being set to 1.0
//  Also, I think an error has been introduced here.
//  Looks like the backbone nodes are always being set to 1.0
//  but this must be checked.
//  Is this correctly handling the initial sampling probabilities
//    for the backbone of the tree? Must check Ei settings.
//
//  DLR comment 12.15.2014
//      Not sure about above comment, looks OK to me.
//

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
            
            
            double LL = computeSpExProbBranch(node->getLfDesc());
            double LR = computeSpExProbBranch(node->getRtDesc());
            
            //std::cout << "main: " << LL << "\t" << LR << std::endl;
            
            logLikelihood += LL + LR;
            
            //logLikelihood += computeSpExProbBranch(node->getLfDesc());
            
 
            //logLikelihood += computeSpExProbBranch(node->getRtDesc());
 
            // Does not include root node, so it is conditioned
            // on basal speciation event occurring:
            if (node != _tree->getRoot()) {
                logLikelihood  += log(node->getNodeLambda());
                node->setDinit(1.0);
            }
        }
    }

    // COMMENT
    //std::cout << "logLik final pre-pres\t" << logLikelihood << std::endl;

    
    logLikelihood += computePreservationLogProb();
 
    // COMMENT
    //std::cout << "logLik final post-pres\t" << logLikelihood << std::endl;
    
    
    return logLikelihood;
}


double SpExModel::computeSpExProbBranch(Node* node)
{
    // COMMENT
    //std::cout << "NOde ID\t" << node << "isInternal: " << node->isInternal();
    //std::cout << "Dinit: " << node->getDinit() << std::endl;
    
    double logLikelihood = 0.0;

    double D0 = node->getDinit();    // Initial speciation probability
    double E0 = node->getEinit();    // Initial extinction probability
    
    // 3 scenarios:
    //   i. node is extant tip
    //   ii. node is fossil last occurrence, an unsampled or extinct tip
    //   iii. node is internal. Will now treat separately.
    //  case i and iii can be treated the same
    
    
    // Test if node is extant
    bool isExtant = (std::abs(node->getTime() - _observationTime)) < 0.00001;
    
    // COMMENT
    //std::cout << node->getTime() << "\t" << isExtant << std::endl;
    
    /****************/

    
    if (node->isInternal() == false & isExtant == false){
    // case 1: node is fossil tip
    
        double ddt = _observationTime - node->getTime();
        
        double startTime = node->getBrlen() + ddt;
        double endTime = startTime;
        
        while (startTime > node->getBrlen()){
            startTime -= _segLength;
            if (startTime < node->getBrlen() ){
                startTime = node->getBrlen();
            }
            double deltaT = endTime - startTime;
            
            double curLam = node->computeSpeciationRateIntervalRelativeTime
            (startTime, endTime);
            
            double curMu = node->computeExtinctionRateIntervalRelativeTime
            (startTime, endTime);
            
            // TODO: curPsi can be computed once we have interval-specific psi values
            double curPsi = _preservationRate;
            
            double spProb = 0.0;
            double exProb = 0.0;
        
            computeSpExProb(spProb, exProb, curLam, curMu, curPsi, D0, E0, deltaT);
            
            E0 = exProb;
            
            endTime = startTime;
            
        }
        
        // Prob that lineage went extinct before present
        // E0 could be the new D0 for the next calculation
        //  however, we will factor this out and start with 1.0.
        
        // COMMENT
        //std::cout << "exprob\t" << E0 << std::endl;
        
        logLikelihood += log(E0);
        node->setDinit(1.0);
        D0 = 1.0;
        // current value of E0 can now be passed on for
        // calculation down remainder of branch (eg, the observed segement)
 
    }
    
    
    double startTime = node->getBrlen();
    double endTime = node->getBrlen();

    while (startTime > 0) {
        startTime -= _segLength;
        if (startTime < 0) {
            startTime = 0.0;
        }

        // Important param for E0 calculations:
        int n_events = node->getBranchHistory()->getNumberOfEventsOnInterval(startTime, endTime);
        
        double deltaT = endTime - startTime;

        double curLam = node->computeSpeciationRateIntervalRelativeTime
            (startTime, endTime);

        double curMu = node->computeExtinctionRateIntervalRelativeTime
            (startTime, endTime);

        double curPsi = _preservationRate;

        
        double spProb = 0.0;
        double exProb = 0.0;

        // Compute speciation and extinction probabilities and store them
        // in spProb and exProb (through reference passing)
        computeSpExProb(spProb, exProb, curLam, curMu, curPsi, D0, E0, deltaT);

        logLikelihood += std::log(spProb);
        
        // COMMENT
        //std::cout << node << "\t" << spProb << std::endl;
        
        D0 = 1.0;
        
        // E0 calculations:
        //     See if branch event occurred on segment.
        //     If no event occurred:
        //         If still on branch
        //             set next E0 to exProb
        //         If reached end of branch:
        //              Set E0 for parent node to exProb
        //     If event occurred on segment:
        //         Recompute E0 from next tstart, using
        //           node event from ancestor
        //         If reached end of branch:
        //           set ancestral nodel to this new E0
        
        if (n_events == 0){
            E0 = exProb;
        }else{
            // recompute E0.
            // A bit of a pain as we have to
            //   redo this from the present backwards in time, piecewise.
            // Should be same calculation as for case above where node
            //   is a fossil tip.
            
            double ddt = _observationTime - node->getTime();
            double st = node->getBrlen() + ddt;
            double et = st;
            
            E0 = node->getEtip();
            
            while (st > startTime){
                
                st -= _segLength;
                if (st <= startTime){
                    st = startTime;
                }
                double tt = et - st;
                
                double clam = node->computeSpeciationRateIntervalRelativeTime
                (st, et);
                
                double cmu = node->computeExtinctionRateIntervalRelativeTime
                (st, et);
                
                double cpsi = _preservationRate;
                
                double sprob = 0.0;
                double eprob = 0.0;
    
                computeSpExProb(sprob, eprob, clam, cmu, cpsi, (double)1.0, E0, tt);
                
                E0 = eprob;
                et = st;
        
            }
            
            
        }
    
        // E0 = exProb;

        endTime = startTime;
    }
    
    // What to use as E0 for the start of next downstream branch?
    // At this point, E0 is actually the E0 for the parent node.
    // Here, we will arbitrarily take this to be the left extinction value
    // Will be identical for right and left nodes at this point
    //      e.g., E0 at parent node coming from right or left descendant

    Node* parent = node->getAnc();
    if (node == parent->getLfDesc()) {
        parent->setEinit(E0);
    }

    // COMMENT
    //std::cout << node->getBrlen() << "\t" << logLikelihood << std::endl;
    //std::cout << "ExProb: " << E0 << "\tLogLik\t" << logLikelihood << std::endl;
    
    
    // Clearly a problem if extinction values approaching/equaling 1.0
    // If so, set to -Inf, leading to automatic rejection of state
    if (E0 > _extinctionProbMax) {
        return -INFINITY;
    }

    // To CONDITION on survival, uncomment the line below;
    // or to NOT condition on survival, comment the line below
    if (parent == _tree->getRoot()) {
        
        // TODO: E0 should probably be E0^2  in this equation
        
        logLikelihood -= std::log(1.0 - E0);
    }
    
    return logLikelihood;
}


//  Notes for the fossil process:
//
//      
//
//

////////////  Notes for the non-fossil process
//
//  If we let:
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
//

// TODO:: Full document equations.
//          equations above are for reconstructed process only.
//          equations as implemented allow for fossilized process.

void SpExModel::computeSpExProb(double& spProb, double& exProb,
                                double lambda, double mu, double psi, double D0, double E0, double deltaT)
{
    
    double FF = lambda - mu - psi;
    double c1 = std::abs(std::sqrt( FF * FF  + (4.0 * lambda * psi) ));
    double c2 = (-1.0) * (FF - 2.0 * lambda * (1.0 - E0)) / c1;
    
    double A = std::exp((-1.0) * c1 * deltaT) * (1.0 - c2);
    double B = c1 * (A - (1 + c2)) / (A + (1.0 + c2));
    
    exProb = (lambda + mu + psi + B)/ (2.0 * lambda);
    
    // splitting up the speciation calculation denominator:
    
    double X = std::exp(c1 * deltaT) * (1.0 + c2)*(1.0 + c2);
    double Y = std::exp((-1.0) * c1 * deltaT) * (1.0 - c2) * (1.0 - c2);
    
    spProb = (4.0 * D0) / ( (2.0 * ( 1 - (c2 * c2)) ) + X + Y );

}

/*
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
*/







double SpExModel::computeLogPrior()
{
    double logPrior = 0.0;

    SpExBranchEvent* rootEvent = static_cast<SpExBranchEvent*>(_rootEvent);

    logPrior += _prior.lambdaInitRootPrior(rootEvent->getLamInit());
    if (rootEvent->isTimeVariable()) {
        logPrior += _prior.lambdaShiftRootPrior(rootEvent->getLamShift());
    }
    logPrior += _prior.muInitRootPrior(rootEvent->getMuInit());
    logPrior += _prior.muShiftRootPrior(rootEvent->getMuShift());

    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        SpExBranchEvent* event = static_cast<SpExBranchEvent*>(*it);

        logPrior += _prior.lambdaInitPrior(event->getLamInit());
        if (event->isTimeVariable()) {
            logPrior += _prior.lambdaShiftPrior(event->getLamShift());
        }
        logPrior += _prior.muInitPrior(event->getMuInit());
        logPrior += _prior.muShiftPrior(event->getMuShift());
    }

    // Here's prior density on the event rate
    logPrior += _prior.poissonRatePrior(_eventRate);

    return logPrior;
}


double SpExModel::computePreservationLogProb()
{
    double logLik = _numberOccurrences * std::log(_preservationRate);
    return logLik;
}




