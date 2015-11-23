#include "Model.h"
#include "Random.h"
#include "Settings.h"
#include "Prior.h"
#include "Tree.h"
#include "EventNumberProposal.h"
#include "EventNumberForBranchProposal.h"
#include "MoveEventProposal.h"
#include "EventRateProposal.h"
 

#include "Node.h"
#include "BranchEvent.h"
#include "BranchHistory.h"
#include "Tools.h"

#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>

#define ENABLE_HASTINGS_RATIO_BUG


Model::Model(Random& random, Settings& settings) :
    _random(random), _settings(settings), _prior(_random, &_settings),
    _tree(new Tree(_random, _settings))
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

    // Add proposals
    _proposals.push_back(new EventNumberProposal(random, settings, *this));
    _proposals.push_back
        (new EventNumberForBranchProposal(random, settings, *this));
    _proposals.push_back(new MoveEventProposal(random, settings, *this));
    _proposals.push_back
        (new EventRateProposal(random, settings, *this, _prior));

}


Model::~Model()
{
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        delete *it;
    }

    delete _tree;

    // Delete all proposals (including those created by derived classes)
    for (Proposal* proposal : _proposals) {
        delete proposal;
    }
    
    delete _rootEvent;
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
            setMeanBranchParameters();
        } else {
            double deltaT = x->getTime() - eventTime;
            double newMapTime = x->getMapStart() + deltaT;

            BranchEvent* newEvent =
                newBranchEventWithReadParameters(x, newMapTime, parameters);
            addEventToTree(newEvent);
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
    // Add un-normalized weights of proposals
    for (Proposal* proposal : _proposals) {
        _updateWeights.push_back(proposal->weight());
    }

    // Sum all weights
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


void Model::proposeNewState()
{
    int parameterToUpdate = chooseParameterToUpdate();
    _lastParameterUpdated = parameterToUpdate;

    Proposal* proposal = _proposals[parameterToUpdate];
    proposal->propose();

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

BranchEvent* Model::addFixedParameterEventToRandomLocation()
{
    double aa = _tree->getRoot()->getMapStart();
    double bb = _tree->getTotalMapLength();
    double x = _random.uniform(aa, bb);
    
    BranchEvent* newEvent = newBranchEventWithParametersFromSettings(x);
    
    return addEventToTree(newEvent);
}


BranchEvent* Model::addRandomEventToTreeOnRandomBranch()
{
    Node* randomNode = _tree->getRandomNonRootNode();

    double mapStart = randomNode->getMapStart();
    double mapEnd = randomNode->getMapEnd();

    double randomMapPoint = _random.uniform(mapStart, mapEnd);
    BranchEvent* newEvent = newBranchEventWithRandomParameters(randomMapPoint);

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


// Model::isEventConfigurationValid
//   Tests whether an event can be placed in a triad configuration:
//   This is when you have an event on a particular branch and on both
//   descendant branches. The configuration leads to an unusual pathology whereby
//   speciation rate can load up on the node bracketed by events.
//   This logical check can be used to reject these configurations
//   as they arise if desired.
//   Define this configuration as having ancestral branch XX
//   and YY and ZZ denoting right and left descendant branches, respectively

//   TODO: fix validateEventConfiguration
//   Oct 2015: First of all, this is technically inappropriate
//              but it is still an important option.
//              So, from one perspective, it is a "Hastings Ratio Bug".
//              However, I think the code itself is not doing what it is supposed
//              to be doing.
//              It should (if iplemented correctly) have a theoretically very small effect
//              on the posterior. It has quite a large effect, suggesting it is inappropriately
//              rejecting many other event configurations besides the "triad" configuration
//              described above.

bool Model::isEventConfigurationValid(BranchEvent* be)
{
    
#ifndef ENABLE_HASTINGS_RATIO_BUG
    std::cout << "\n********* ERROR *********** " << std::endl;
    std::cout << "validateEventConfiguration is no longer a valid option in BAMM" << std::endl;
    std::cout << "\nYou will have to go to src/Model.cpp and recompile " << std::endl;
    std::cout << "after uncommenting the ENABLE_HASTINGS_RATIO_BUG macro" << std::endl;
    std::cout << "\n\n" << std::endl;
    
    exit(0);
    
#endif
    
    bool isValidConfig = false;

    bool forwardConfigValid = false;
    bool backwardConfigValid = false;
    
    if (be->getEventNode() == _tree->getRoot()) {
        Node* rt = _tree->getRoot()->getRtDesc();
        Node* lf = _tree->getRoot()->getLfDesc();
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0 &&
            lf->getBranchHistory()->getNumberOfBranchEvents() > 0) {
            // Events on both descendants of root. This fails.
            isValidConfig = false;
        } else {
            isValidConfig = true;
        }

    } else {
        int badsum = 0;

        // Backward check phase: tests whether events on YY or ZZ are valid
        
        
        Node* anc = be->getEventNode()->getAnc();
        Node* lf = anc->getLfDesc();
        Node* rt = anc->getRtDesc();

        
    
        // Test ancestor for events on branch

        if (anc == _tree->getRoot()) {
            badsum++;
        } else if (anc->getBranchHistory()->getNumberOfBranchEvents() > 0) {
            badsum++;
        } else {
            // nothing
        }

        // Test lf desc
        if (lf->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        // Test rt desc
        if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0)
            badsum++;

        if (badsum == 3) {
            backwardConfigValid = false;
        } else if (badsum < 3) {
            backwardConfigValid = true;
        } else {
            log(Error) << "Problem in Model::isEventConfigurationValid\n";
            std::exit(1);
        }
        
        // Forward check phase: tests whether events on XX are valid
        // does not count focal branch since there is
        // obviously an event on branch defined by event be
        anc = be->getEventNode();
        lf = anc->getLfDesc();
        rt = anc->getRtDesc();
        badsum = 0;
        
        if (lf != NULL && rt != NULL && backwardConfigValid){
            if (lf->getBranchHistory()->getNumberOfBranchEvents() > 0){
                badsum++;
            }
            if (rt->getBranchHistory()->getNumberOfBranchEvents() > 0){
                badsum++;
            }
            
            if (badsum == 2) {
                forwardConfigValid = false;
            } else{
                forwardConfigValid = true;
            }
    
        }
        
        if (forwardConfigValid && backwardConfigValid){
            isValidConfig = true;
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

void Model::printEventValidStatus()
{
    
    std::cout << "\nChecking event configuration: Model::printEventValidStatus " << std::endl;
    
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it) {
        
        bool isValid = isEventConfigurationValid((*it));

        std::cout << (*it)->getEventNode() << "\tisValid:\t " << isValid << std::endl;
    }


}



bool Model::testEventConfigurationComprehensive()
{
    bool isValidAll = true;
    
    EventSet::iterator it;
    for (it = _eventCollection.begin(); it != _eventCollection.end(); ++it){
        bool isValidSingle = isEventConfigurationValid((*it));
        if (!isValidSingle){
            isValidAll = false;
            break;
        }
    }
    return isValidAll;
}













