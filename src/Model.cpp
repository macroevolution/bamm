#include "Model.h"
#include "Random.h"
#include "Tree.h"
#include "Settings.h"
#include "Prior.h"
#include "Node.h"
#include "BranchEvent.h"
#include "BranchHistory.h"
#include "Log.h"
#include "Tools.h"

#include <string>
#include <fstream>
#include <cstdlib>


Model::Model(Random& random, Settings& settings) :
    _random(random), _settings(settings), _prior(_random, &_settings),
    _tree(new Tree(_random, _settings)),
    _eventNumberProposal(random, settings, *this),
    _moveEventProposal(random, settings, *this),
    _eventRateProposal(random, settings, *this, _prior)
{
    // Initialize event rate to generate expected number of prior events
    _eventRate = 1 / _settings.get<double>("poissonRatePrior");

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
}


Model::~Model()
{
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        delete *it;
    }

    delete _tree;
}


void Model::initializeModelFromEventDataFile(const std::string& fileName)
{
    std::ifstream inputFile(fileName.c_str());

    if (!inputFile) {
        log(Error) << "Could not read event data file "
            << "<<" << fileName << ">>.\n";
        std::exit(1);
    }

    log() << "Initializing model from <<" << fileName << ">>...\n";

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(inputFile, line)) {
        lines.push_back(line);
    }

    inputFile.close();

    int eventCount = 0;
    int prevGeneration = 0;

    // The relevant events are at the bottom of the file (last generation)
    for (int i = lines.size() - 1; i != -1; --i) {
        const std::vector<std::string>& tokens = split_string(lines[i], ',');

        // Get the generation, but if it differs from previous, stop
        int gen = convert_string<int>(tokens[0]);
        if (prevGeneration == 0) {
            prevGeneration = gen;
        } else {
            if (gen != prevGeneration) {
                break;
            }
        }

        std::string species_1 = tokens[1];
        std::string species_2 = tokens[2];
        double eventTime = convert_string<double>(tokens[3]);

        std::vector<std::string> parameters;
        for (int i = 4; i < (int)tokens.size(); ++i) {
            parameters.push_back(tokens[i]);
        }

        Node* x = NULL;
        
        if ((species_1 != "NA") && (species_2 != "NA")) {
            x = _tree->getNodeMRCA(species_1, species_2);
        } else if ((species_1 != "NA") && (species_2 == "NA")) {
            x = _tree->getNodeByName(species_1);
        } else {
            log(Error) << "Either both species are NA or the second species "
                << "is NA\nwhile reading the event data file.";
            std::exit(1);
        }

        if (x == _tree->getRoot()) {
            // Set the root event with model-specific parameters
            setRootEventWithReadParameters(parameters);
        } else {
            double deltaT = x->getTime() - eventTime;
            double newMapTime = x->getMapStart() + deltaT;

            BranchEvent* newEvent =
                newBranchEventWithReadParameters(x, newMapTime, parameters);
            newEvent->getEventNode()->getBranchHistory()->
                addEventToBranchHistory(newEvent);

            _eventCollection.insert(newEvent);
            forwardSetBranchHistories(newEvent);
            setMeanBranchParameters();
        }

        eventCount++;
    }

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
    _updateWeights.push_back(_settings.get<double>("updateRateEventNumber"));
    _updateWeights.push_back(_settings.get<double>("updateRateEventPosition"));
    _updateWeights.push_back(_settings.get<double>("updateRateEventRate"));

    // Defined by derived class
    initializeSpecificUpdateWeights();
}


void Model::proposeNewState()
{
    int parameterToUpdate = chooseParameterToUpdate();
    _lastParameterUpdated = parameterToUpdate;

    Proposal* proposal = NULL;

    if (parameterToUpdate == 0) {
        proposal = &_eventNumberProposal;
    } else if (parameterToUpdate == 1) {
        proposal = &_moveEventProposal;
    } else if (parameterToUpdate == 2) {
        proposal = &_eventRateProposal;
    } else {
        // Defined by derived class
        proposal = getSpecificProposal(parameterToUpdate);
    }

    if (proposal != NULL) {
        proposal->propose();
    }

    _lastProposal = proposal;
}


int Model::chooseParameterToUpdate()
{
    double r = _random.uniform();

    for (int i = 0; i < (int)_updateWeights.size(); i++) {
        if (r < _updateWeights[i]) {
            return i;
        }
    }

    return -1;
}


BranchEvent* Model::addRandomEventToTree()
{
    double aa = _tree->getRoot()->getMapStart();
    double bb = _tree->getTotalMapLength();
    double x = _random.uniform(aa, bb);
                
    BranchEvent* newEvent = newBranchEventWithRandomParameters(x);
    return addEventToTree(newEvent);
}


BranchEvent* Model::addEventToTree(BranchEvent* newEvent)
{
    // Add the event to the branch history.
    // Always done after event is added to tree.
    newEvent->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(newEvent);
                
    _eventCollection.insert(newEvent);
    forwardSetBranchHistories(newEvent);
    setMeanBranchParameters();
                            
    _lastEventModified = newEvent;

    return newEvent;
}


BranchEvent* Model::chooseEventAtRandom(bool includeRoot)
{   
    EventSet& events = _eventCollection;
    int numberOfEvents = (int)events.size();

    int eventIndex = 0;
    if (includeRoot) {
        eventIndex = _random.uniformInteger(0, numberOfEvents);
    } else {
        eventIndex = _random.uniformInteger(0, numberOfEvents - 1);
    }
    
    if (eventIndex == numberOfEvents) {
        return getRootEvent();
    } else {
        EventSet::iterator it = events.begin();
        std::advance(it, eventIndex);
        return *it;
    }
}


BranchEvent* Model::removeEventFromTree(BranchEvent* be)
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

    forwardSetBranchHistories(newLastEvent);

    setMeanBranchParameters();

    return be;
}


BranchEvent* Model::removeRandomEventFromTree()
{
    int numEvents = (int)_eventCollection.size();

    // Can only delete event if more than root node present.
    if (numEvents == 0) {
        return NULL;
    }

    int counter = 0;
    double xx = _random.uniform();
    int chosen = (int)(xx * (double)numEvents);

    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        if (counter++ == chosen) {
            return removeEventFromTree(*it);
        }
    }

    return NULL;
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


void Model::acceptProposal()
{
    if (_lastProposal != NULL) {
        _lastProposal->accept();
        _acceptCount++;
        _acceptLast = 1;
    } else {
        _acceptLast = -1;
    }
}


void Model::rejectProposal()
{
    if (_lastProposal != NULL) {
        _lastProposal->reject();
        _rejectCount++;
        _acceptLast = 0;
    } else {
        _acceptLast = -1;
    }
}


double Model::acceptanceRatio()
{
    if (_lastProposal == NULL) {
        return 0.0;
    }

    return _lastProposal->acceptanceRatio();
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


double Model::safeExponentiation(double x)
{
    if (x > 0.0)
        return 1.0;
    else if (x < -100.0)
        return 0.0;
    else
        return std::exp(x);
}


void Model::setTemperatureMH(double x)
{
    if ((x >= 0) && (x <= 1.0)) {
        _temperatureMH = x;
    } else {
        log(Error) << "Attempt to set invalid temperature in "
            << "Model::setModelTemperature: " << x << "\n";
        std::exit(1);
    }
}
