#include "Model.h"
#include "MbRandom.h"
#include "Tree.h"
#include "Settings.h"
#include "Prior.h"
#include "Node.h"
#include "BranchEvent.h"
#include "BranchHistory.h"
#include "Log.h"

#include <string>
#include <fstream>
#include <cstdlib>


Model::Model(MbRandom* rng, Tree* tree, Settings* settings, Prior* prior) :
    _rng(rng), _tree(tree), _settings(settings), _prior(prior)
{
    // Reduce weird autocorrelation of values at start by calling RNG
    // a few times. TODO: Why is there a weird autocorrelation?
    for (int i = 0; i < 100; i++)
        _rng->uniformRv();

    // Event location scale is relative to the maximum root-to-tip length
    _scale = _settings->getUpdateEventLocationScale() *
        _tree->maxRootToTipLength();

    _updateEventRateScale = _settings->getUpdateEventRateScale();
    _localGlobalMoveRatio = _settings->getLocalGlobalMoveRatio();
    
    _poissonRatePrior = _settings->getPoissonRatePrior();

    // Initialize event rate to generate expected number of prior events
    _eventRate = 1 / _settings->getPoissonRatePrior();

    _acceptCount = 0;
    _rejectCount = 0;
    _acceptLast = -1;

    _lastDeletedEventMapTime = 0;

    _logQRatioJump = 0.0;

    // Initial setting for temperature = 1.0
    // This must be explicitly set by calling the public
    // function Model::setModelTemperature
     
    _temperatureMH = 1.0;
}


// This method needs to be called by the derived class
void Model::finishConstruction()
{
    calculateUpdateWeights();

    int initialNumberOfEvents = _settings->getInitialNumberEvents();
    for (int i = 0; i < initialNumberOfEvents; i++) {
        addEventToTree();
    }
}


Model::~Model()
{
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        delete *it;
    }
}


void Model::initializeModelFromEventDataFile()
{
    std::string inputFileName(_settings->getEventDataInfile());
    std::ifstream inputFile(inputFileName.c_str());

    if (!inputFile.good()) {
        log(Error) << "<<" << inputFileName << ">> is a bad file name.\n";
        std::exit(1);
    }

    log() << "Initializing model from <<" << inputFileName << ">>\n";

    std::string species1;
    std::string species2;
    double eTime;

    int eventCount = 0;
    while (inputFile) {
        inputFile >> species1;
        inputFile >> species2;
        inputFile >> eTime;

        // Read the model-specific parameters
        readModelSpecificParameters(inputFile);

        // TODO: Might need to getline here to read last \n

        Node* x = NULL;
        
        if ((species1 != "NA") && (species2 != "NA")) {
            x = _tree->getNodeMRCA(species1.c_str(), species2.c_str());
        } else if ((species1 != "NA") && (species2 == "NA")) {
            x = _tree->getNodeByName(species1.c_str());
        } else {
            log(Error) << "Either both species are NA or the second species "
                << "is NA\nwhile reading the event data file.";
            std::exit(1);
        }

        if (x == _tree->getRoot()) {
            // Set the root event with model-specific parameters
            setRootEventWithReadParameters();
        } else {
            double deltaT = x->getTime() - eTime;
            double newMapTime = x->getMapStart() + deltaT;

            BranchEvent* newEvent =
                newBranchEventWithReadParameters(x, newMapTime);
            newEvent->getEventNode()->getBranchHistory()->
                addEventToBranchHistory(newEvent);

            _eventCollection.insert(newEvent);
            forwardSetBranchHistories(newEvent);
            setMeanBranchParameters();
        }

        eventCount++;
    }

    inputFile.close();

    log() << "Read a total of " << eventCount << " events.\n";
    log() << "Added " << _eventCollection.size() << " "
          << "pre-defined events to tree, plus root event.\n";
}


void Model::forwardSetBranchHistories(BranchEvent* x)
{
    // If there is another event occurring more recent (closer to tips),
    // do nothing. Even just sits in BranchHistory but doesn't affect
    // state of any other nodes.

    // This seems circular, but what else to do?
    // given an event (which references the node defining the branch on which
    // event occurs) you get the corresponding branch history and the last
    // event since the events will have been inserted in the correct order.

    Node* myNode = x->getEventNode();

    if (x == _rootEvent) {
        forwardSetHistoriesRecursive(myNode->getLfDesc());
        forwardSetHistoriesRecursive(myNode->getRtDesc());
    } else if (x == myNode->getBranchHistory()->getLastEvent()) {
        // If true, x is the most tip-wise event on branch.
        myNode->getBranchHistory()->setNodeEvent(x);

        // If myNode is not a tip
        if (myNode->getLfDesc() != NULL && myNode->getRtDesc() != NULL) {
            forwardSetHistoriesRecursive(myNode->getLfDesc());
            forwardSetHistoriesRecursive(myNode->getRtDesc());
        }
        // Else: node is a tip; do nothing
    }
    // Else: there is another more tipwise event on the same branch; do nothing
}


/*
    If this works correctly, this will take care of the following:
    1. if a new event is created or added to tree,
       this will forward set all branch histories from the insertion point
    2. If an event is deleted, you find the next event rootwards,
       and call forwardSetBranchHistories from that point. It will replace
       settings due to the deleted node with the next rootwards node.
*/

void Model::forwardSetHistoriesRecursive(Node* p)
{
    // Get event that characterizes parent node
    BranchEvent* lastEvent = p->getAnc()->getBranchHistory()->getNodeEvent();

    // Set the ancestor equal to the event state of parent node:
    p->getBranchHistory()->setAncestralNodeEvent(lastEvent);

    // Ff no events on the branch, go down to descendants and do same thing;
    // otherwise, process terminates (because it hits another event on branch
    if (p->getBranchHistory()->getNumberOfBranchEvents() == 0) {
        p->getBranchHistory()->setNodeEvent(lastEvent);

        if (p->getLfDesc() != NULL) {
            forwardSetHistoriesRecursive(p->getLfDesc());
        }

        if (p->getRtDesc() != NULL) {
            forwardSetHistoriesRecursive(p->getRtDesc());
        }
    }
}


void Model::calculateUpdateWeights()
{
    initializeUpdateWeights();

    // Add all weights
    double sumWeights = 0.0;
    for (int i = 0; i < (int)_updateWeights.size(); i++) {
        sumWeights += _updateWeights[i];
    }

    // Calculate cumulative weights
    for (int i = 1; i < (int)_updateWeights.size(); i++) {
        _updateWeights[i] += _updateWeights[i - 1];
    }

    // Normalize by the sum of all weights
    for (int i = 0; i < (int)_updateWeights.size(); i++) {
        _updateWeights[i] /= sumWeights;
    }
}


void Model::initializeUpdateWeights()
{
    _updateWeights.push_back(_settings->getUpdateRateEventNumber());
    _updateWeights.push_back(_settings->getUpdateRateEventPosition());
    _updateWeights.push_back(_settings->getUpdateRateEventRate());

    // Defined by derived class
    initializeSpecificUpdateWeights();
}


void Model::proposeNewState()
{
    int parameterToUpdate = chooseParameterToUpdate();
    _lastParameterUpdated = parameterToUpdate;

    if (parameterToUpdate == 0) {
        changeNumberOfEventsMH();
    } else if (parameterToUpdate == 1) {
        moveEventMH();
    } else if (parameterToUpdate == 2) {
        updateEventRateMH();
    } else {
        // Defined by derived class
        proposeSpecificNewState(parameterToUpdate);
    }
}


int Model::chooseParameterToUpdate()
{
    double r = _rng->uniformRv();

    for (int i = 0; i < (int)_updateWeights.size(); i++) {
        if (r < _updateWeights[i]) {
            return i;
        }
    }

    return -1;
}


void Model::addEventToTree()
{
    double aa = _tree->getRoot()->getMapStart();
    double bb = _tree->getTotalMapLength();
    double x = _rng->uniformRv(aa, bb);
                
    addEventToTree(x);
}


// Adds event to tree based on reference map value
// - Adds to branch history set
// - Inserts into _eventCollection

void Model::addEventToTree(double x)
{
    BranchEvent* newEvent = newBranchEventWithRandomParameters(x);
            
    // Add the event to the branch history.
    // Always done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(newEvent);
                
    _eventCollection.insert(newEvent);
    forwardSetBranchHistories(newEvent);
    setMeanBranchParameters();
                            
    _lastEventModified = newEvent;
}


// TODO: ctr is not doing anything
BranchEvent* Model::chooseEventAtRandom()
{
    int numEvents = (int)_eventCollection.size();

    if (numEvents == 0) {
        log(Error) << "Number of events is zero.\n";
        std::exit(1);
    }

    int ctr = 0;
    double xx = _rng->uniformRv();
    int chosen = (int)(xx * (double)numEvents);

    EventSet::iterator sit = _eventCollection.begin();

    for (int i = 0; i < chosen; i++) {
        ++sit;
        ctr++;
    }

    return *sit;
}


void Model::eventLocalMove(void)
{
    eventMove(true);
}


void Model::eventGlobalMove(void)
{
    eventMove(false);
}


// If events are on tree: choose event at random,
// move locally (or globally) and forward set branch histories etc.
// Should also store previous event information to revert to previous

// If parameter local == true, does a local move;
// otherwise, it does a global move

void Model::eventMove(bool local)
{
    if (getNumberOfEvents() > 0) {
        // The event to be moved
        BranchEvent* chosenEvent = chooseEventAtRandom();

        // This is the event preceding the chosen event:
        // histories should be set forward from here.
        BranchEvent* previousEvent = chosenEvent->getEventNode()->
            getBranchHistory()->getLastEvent(chosenEvent);

        // Set this history variable in case move is rejected
        _lastEventModified = chosenEvent;

        chosenEvent->getEventNode()->getBranchHistory()->
            popEventOffBranchHistory(chosenEvent);

        if (local) {
            // Get step size for move
            double step = _rng->uniformRv(0, _scale) - 0.5 * _scale;
            chosenEvent->moveEventLocal(step);
        } else {
            chosenEvent->moveEventGlobal();
        }

        chosenEvent->getEventNode()->getBranchHistory()->
            addEventToBranchHistory(chosenEvent);

        // Get last event from the theEventNode, forward set its history.
        // Then go to the "moved" event and forward set its history.

        forwardSetBranchHistories(previousEvent);
        forwardSetBranchHistories(chosenEvent);
    }

    setMeanBranchParameters();
}


// Used to reset position of event if move is rejected

void Model::revertMovedEventToPrevious()
{
    // Get last event from position of event to be removed
    BranchEvent* newLastEvent = _lastEventModified->getEventNode()->
        getBranchHistory()->getLastEvent(_lastEventModified);

    // Pop event off its new position
    _lastEventModified->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_lastEventModified);

    // Reset nodeptr, reset mapTime
    _lastEventModified->revertOldMapPosition();

    // Now reset forward from _lastEventModified (new position)
    // and from newLastEvent, which holds 'last' event before old position
    _lastEventModified->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(_lastEventModified);

    // Forward set from new position
    forwardSetBranchHistories(newLastEvent);

    // Forward set from event immediately rootwards from previous position
    forwardSetBranchHistories(_lastEventModified);

    // Set _lastEventModified to NULL because it has already been reset.
    // Future implementations should check whether this is NULL
    // before attempting to use it to set event

    _lastEventModified = NULL;

    setMeanBranchParameters();
}


void Model::deleteEventFromTree(BranchEvent* be)
{
    if (be ==_rootEvent) {
        log(Error) << "Can't delete root event.\n";
        std::exit(1);
    }

    // Erase from branch history
    Node* currNode = be->getEventNode();

    // Get event downstream of i
    BranchEvent* newLastEvent = currNode->getBranchHistory()->getLastEvent(be);

    _lastDeletedEventMapTime = be->getMapTime();

    setDeletedEventParameters(be);
    _logQRatioJump = calculateLogQRatioJump();

    currNode->getBranchHistory()->popEventOffBranchHistory(be);

    // Cannot remove "be" with _eventCollection.erase(be) because
    // it is not always found in the collection, even though it is there.
    // It seems that, because event pointers are compared by comparing
    // doubles, the event is sometimes not seen (floating-point issue).
    bool eventFound = false;
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        if (*it == be) {    // Compare pointers directly, not using comparer
            _eventCollection.erase(it);
            eventFound = true;
            break;
        }
    }

    if (!eventFound) {
        log(Error) << "Could not find event to delete.\n";
        std::exit(1);
    }

    delete be;

    forwardSetBranchHistories(newLastEvent);

    setMeanBranchParameters();
}


void Model::deleteRandomEventFromTree()
{
    int numEvents = (int)_eventCollection.size();

    // Can only delete event if more than root node present.
    if (numEvents == 0) {
        return;
    }

    int counter = 0;
    double xx = _rng->uniformRv();
    int chosen = (int)(xx * (double)numEvents);

    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        if (counter++ == chosen) {
            deleteEventFromTree(*it);
            break;
        }
    }
}


void Model::restoreLastDeletedEvent()
{   
    // Use constructor for speciation and extinction
    BranchEvent* newEvent = newBranchEventFromLastDeletedEvent();

    // Add the event to the branch history.
    // Always done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(newEvent);

    _eventCollection.insert(newEvent);

    // Event is now inserted into branch history;
    // however, branch histories must be updated.
    forwardSetBranchHistories(newEvent);

    setMeanBranchParameters();
}


void Model::changeNumberOfEventsMH()
{
    bool addEvent = _rng->uniformRv() < 0.5;
    if (addEvent || _eventCollection.size() == 0) {
        addEventMH();
    } else {
        removeEventMH();
    }
}


void Model::addEventMH()
{
    double oldLogLikelihood = getCurrentLogLikelihood();
    double oldLogPrior = computeLogPrior();

    int K = (int)_eventCollection.size();
    double qRatio = (K > 0) ? 1.0 : 0.5;

    addEventToTree();
    setMeanBranchParameters();

    double logLikelihood = computeLogLikelihood();
    double logPrior = computeLogPrior();

    double logHR = computeEventGainLogHR(K, logLikelihood, oldLogLikelihood,
        logPrior, oldLogPrior, qRatio);

    bool acceptMove = false;
    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    if (!std::isinf(logLikelihoodRatio)) {  // TODO: When is this ever inifinte?
        acceptMove = acceptMetropolisHastings(logHR);
    }

    bool isValidConfig = isEventConfigurationValid(_lastEventModified);
    if (acceptMove && isValidConfig) {
        setCurrentLogLikelihood(logLikelihood);
        _acceptCount++;
        _acceptLast = 1;
    } else {
        deleteEventFromTree(_lastEventModified);
        setMeanBranchParameters();
        _rejectCount++;
        _acceptLast = 0;
    }
}


double Model::computeEventGainLogHR(double K, double logLikelihood,
    double oldLogLikelihood, double logPrior, double oldLogPrior, double qRatio)
{
    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    double logPriorRatio = logPrior - oldLogPrior + std::log(_eventRate) - std::log(K + 1.0);

    double logQratio = std::log(qRatio) - _logQRatioJump;
    
    double logHR = computeLogHastingsRatio(logLikelihoodRatio, logPriorRatio, logQratio);
    
    return logHR;
}


void Model::removeEventMH()
{
    double oldLogLikelihood = getCurrentLogLikelihood();
    double oldLogPrior = computeLogPrior();

    int K = (int)_eventCollection.size();
    double qRatio = (K != 1) ? 1.0 : 2.0;

    deleteRandomEventFromTree();
    setMeanBranchParameters();

    double logLikelihood = computeLogLikelihood();
    double logPrior = computeLogPrior();

    // First get prior ratio:
    double logHR = computeEventLossLogHR(K, logLikelihood, oldLogLikelihood,
        logPrior, oldLogPrior, qRatio);

    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    bool acceptMove = false;
    if (!std::isinf(logLikelihoodRatio)) {
        acceptMove = acceptMetropolisHastings(logHR);
    }

    if (acceptMove) {
        setCurrentLogLikelihood(logLikelihood);
        _acceptCount++;
        _acceptLast = 1;
    } else {
        restoreLastDeletedEvent();
        _rejectCount++;
        _acceptLast = 0;
    }
}


double Model::computeEventLossLogHR(double K, double logLikelihood,
    double oldLogLikelihood, double logPrior, double oldLogPrior, double qRatio)
{

    double logLikelihoodRatio = logLikelihood - oldLogLikelihood;
    
    // This is probability of going from K to K - 1
    // So, prior ratio includes (K / eventRate)
    
    double logPriorRatio = logPrior - oldLogPrior + std::log(K) - std::log(_eventRate);
    
    // For our purposes here, qRatio will include both the proposal ratio
    // and the jumping density of the bijection between parameter spaces
    
    double logQratio = std::log(qRatio) + _logQRatioJump;
    
    double logHR = computeLogHastingsRatio(logLikelihoodRatio, logPriorRatio, logQratio);
    
    return logHR;
}


void Model::moveEventMH()
{
    // Consider proposal rejected (can't move nonexistent event)
    if (_eventCollection.size() == 0) {
        _rejectCount++;
        _acceptLast = 0;
        return;
    }

    double localMoveProb = _localGlobalMoveRatio / (1 + _localGlobalMoveRatio);

    bool isLocalMove = _rng->uniformRv() < localMoveProb;

    if (isLocalMove) {
        // Local move, with event drawn at random
        eventLocalMove();
    } else {
        eventGlobalMove();
    }

    setMeanBranchParameters();

    double logLikelihood = computeLogLikelihood();
    double logLikelihoodRatio = logLikelihood - getCurrentLogLikelihood();
    
    double logHR = computeLogHastingsRatio(logLikelihoodRatio, (double)0.0, double(0.0));
    

    bool isValid = isEventConfigurationValid(_lastEventModified);

    bool acceptMove = false;
    if (!std::isinf(logLikelihoodRatio) && isValid) {
        acceptMove = acceptMetropolisHastings(logHR);
    }

    if (acceptMove) {
        setCurrentLogLikelihood(logLikelihood);
        _acceptCount++;
        _acceptLast = 1;
    } else {
        revertMovedEventToPrevious();
        setMeanBranchParameters();
        _rejectCount++;
        _acceptLast = 0;
    }
}


/*
   Metropolis-Hastings step to update Poisson event rate.
   Note that changing this rate does not affect the likelihood,
   so the priors and qratio determine acceptance rate.
*/

void Model::updateEventRateMH()
{
    double oldEventRate = getEventRate();

    double cterm = std::exp(_updateEventRateScale * (_rng->uniformRv() - 0.5));
    setEventRate(cterm * oldEventRate);

    double logPriorRatio = _prior->poissonRatePrior(getEventRate());
    logPriorRatio -= _prior->poissonRatePrior(oldEventRate);

    double logProposalRatio = std::log(cterm);
    
    double logHR = computeLogHastingsRatio((double)0.0, logPriorRatio, logProposalRatio);
    
    bool acceptMove = acceptMetropolisHastings(logHR);

    if (acceptMove) {
        _acceptCount++;
        _acceptLast = 1;
    } else {
        setEventRate(oldEventRate);
        _rejectCount++;
        _acceptLast = 0;
    }
}


// The original code used here as a template from MrBayes is not correct:
// It was raising the entire numerator of the Hastings Ratio
// to the specified temperature, but temperature should only
// be used for the ratio of posterior probabilities
// WAS:     double r = safeExponentiation( Model::MhColdness * lnR );
// with MhColdness equal to chain temperature (0 < MhColdness <= 1)

bool Model::acceptMetropolisHastings(double lnR)
{
    double r = safeExponentiation( lnR );
    return _rng->uniformRv() < r;
}


bool Model::isEventConfigurationValid(BranchEvent* be)
{
    bool isValidConfig = false;

    if (be->getEventNode() == _tree->getRoot()) {
        Node* rt = _tree->getRoot()->getRtDesc();
        Node* lf = _tree->getRoot()->getLfDesc();
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0 &&
            lf->getBranchHistory()->getNumberOfBranchEvents() > 0) {
            // Events on both descendants of root. This fails.
            isValidConfig = false;
        } else
            isValidConfig = true;

    } else {
        int badsum = 0;

        Node* anc = be->getEventNode()->getAnc();
        Node* lf = anc->getLfDesc();
        Node* rt = anc->getRtDesc();

        // Test ancestor for events on branch

        if (anc == _tree->getRoot())
            badsum++;
        else if (anc->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;
        else {
            // nothing
        }

        // Test lf desc
        if (lf->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        // Test rt desc
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        if (badsum == 3)
            isValidConfig = false;
        else if (badsum < 3)
            isValidConfig = true;
        else {
            log(Error) << "Problem in Model::isEventConfigurationValid\n";
            std::exit(1);
        }
    }

    return isValidConfig;
}


double Model::getMHAcceptanceRate()
{
    return (double)_acceptCount / (_acceptCount + _rejectCount);
}


void Model::resetMHAcceptanceParameters()
{
    _acceptCount = 0;
    _rejectCount = 0;
}


//  Write event data to file for all events "on" tree
//  at a given point in the MCMC chain
    
void Model::getEventDataString(std::stringstream& ss, int generation)
{   
    ss << generation << ",";
    
    BranchEvent* be = _rootEvent;
    Node* xl = _tree->getRoot()->getRandomLeftTipNode();
    Node* xr = _tree->getRoot()->getRandomRightTipNode();
    ss << xl->getName() << ","
       << xr->getName() << ","
       << be->getAbsoluteTime() << ",";
    
    // Implemented in derived class
    getSpecificEventDataString(ss, be);

    if (_eventCollection.size() == 0) {
        return;
    }
    
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        ss << "\n" << generation << ",";
        be = *it;
        if (be->getEventNode()->getLfDesc() == NULL) {
            ss << be->getEventNode()->getName() << "," << "NA" << ",";
        } else {
            Node* xl = be->getEventNode()->getRandomLeftTipNode();
            Node* xr = be->getEventNode()->getRandomRightTipNode();
            ss << xl->getName() << "," << xr->getName() << ",";
        }
        ss << be->getAbsoluteTime() << ",";

        // Implemented in derived class
        getSpecificEventDataString(ss, be);
    }
}


double Model::safeExponentiation(double x)
{
    if (x > 0.0)
        return 1.0;
    else if (x < -100.0)
        return 0.0;
    else
        return std::exp(x);
}


void Model::setModelTemperatureMH(double x)
{
    if ((x >= 0) && (x <= 1.0)) {
        _temperatureMH = x;
    } else {
        log(Error) << "Attempt to set invalid temperature in "
            << "Model::setModelTemperature: " << x << "\n";
        std::exit(1);
    }
}
