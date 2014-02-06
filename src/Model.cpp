#include "Model.h"
#include "MbRandom.h"
#include "Tree.h"
#include "Settings.h"
#include "Prior.h"


double Model::mhColdness = 1.0;


Model::Model(MbRandom* rng, Tree* tree, Settings* settings, Prior* prior) :
    _rng(rng), _tree(tree), _settings(settings), _prior(prior), _gen(0)
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
}


Model::~Model()
{
}

/*
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
*/
