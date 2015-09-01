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
#include <fstream>
   

#define JUMP_VARIANCE_NORMAL 0.05


SpExModel::SpExModel(Random& random, Settings& settings) :
    Model(random, settings)
{
    // Initial values
    _lambdaInit0 = _settings.get<double>("lambdaInit0");
    _lambdaShift0 = _settings.get<double>("lambdaShift0");
    _muInit0 = _settings.get<double>("muInit0");
    _muShift0 = _settings.get<double>("muShift0");

    
    
    // CHECK for paleontological data
    _hasPaleoData = false;
    
    initializeHasPaleoData();
    
    int cs = _settings.get<int>("conditionOnSurvival");
    if (cs == -1){
        if (_hasPaleoData){
            _conditionOnSurvival = false;
        }else{
            _conditionOnSurvival = true;
        }
    }else if (cs == 0){
        _conditionOnSurvival = false;
    }else if (cs == 1){
        _conditionOnSurvival = true;
    }else{
        exitWithError("Invalid initial value for parameter <<conditionOnSurvivial>>");
    }

    // Initialize fossil preservation rate:
    //      will not be relevant if this is not paleo data.
    _preservationRate = _settings.get<double>("preservationRateInit");
    
     
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

    
    //// Change from BranchEvent to SpExBranchEvent:
    BranchEvent* x =  new SpExBranchEvent(_lambdaInit0, _lambdaShift0,
        _muInit0, _muShift0, _initialLambdaIsTimeVariable,
        _tree->getRoot(), _tree, _random, 0);
    
    
    _rootEvent = x;
    _lastEventModified = x;
    
    // DEBUG section
    //SpExBranchEvent* nbe =  static_cast<SpExBranchEvent*>(x);
    //std::cout << "Root event parameters: " << std::endl;
    //std::cout << "\tLambda_init:    " <<  static_cast<SpExBranchEvent*>(x)->getLamInit() << std::endl;
    //std::cout << "\tLambda_shift:   " << static_cast<SpExBranchEvent*>(x)->getLamShift() << std::endl;
    //std::cout << "\tMu_init:        " << static_cast<SpExBranchEvent*>(x)->getMuInit() << std::endl;
    // std::cout << "\tPres_rate:      " << _preservationRate << std::endl;
    
    // DEBUG
    //initializeDebugVectors();
    

    // Set NodeEvent of root node equal to the_rootEvent:
    _tree->getRoot()->getBranchHistory()->setNodeEvent(_rootEvent);

    // Initialize all branch histories to equal the root event
    forwardSetBranchHistories(_rootEvent);

    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();
    
    _extinctionProbMax = _settings.get<double>("extinctionProbMax");

    
    // Initialize by previous event histories (or from initial event number)
    if (_settings.get<bool>("loadEventData")) {
        initializeModelFromEventDataFile(_settings.get("eventDataInfile"));
    } else {
        int initialNumberOfEvents = _settings.get<int>("initialNumberEvents");
        for (int i = 0; i < initialNumberOfEvents; i++) {
            
            // TODO: this adds event to tree with parameters sampled from the prior
            //       should give option to pull parameters from root event (lambdaInit)
            //       etc such that you can add many events but the initial likelihood is
            //       still the same as with a single event
            //addRandomEventToTree();
        
            // This (below) function adds parameters to random locations
            //  but fixes parameters to the initial values specified in control file
            addFixedParameterEventToRandomLocation();
            
            
            // DEBUG
            std::cout << "Event: \t" << i << "\t" << computeLogLikelihood() << std::endl;
            
        }
    }

    if (_settings.get<bool>("validateEventConfiguration")){
        bool isValid = testEventConfigurationComprehensive();
        if (!isValid){
            std::cout << "\nInitial event configuration is invalid\n";
            std::cout << "Please check your output carefully" << std::endl;
            printEventValidStatus();
        }
    }
    
    // DEBUG
    //printEventValidStatus();
    

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

    if (_hasPaleoData){
        // Cannot set this parameter unless you have paleo data....
        _proposals.push_back(new PreservationRateProposal(random, settings, *this, _prior));
    }


    Model::calculateUpdateWeights();
    
    //testPreservationRate();

}


// Evalates settings and tree to determine if
// it is a valid instance of a tree with some paleontological data.
// Sets the _hasPaleoData parameter.
void SpExModel::initializeHasPaleoData()
{
    _numberOccurrences = _settings.get<int>("numberOccurrences");
    
    if (_numberOccurrences > 0 & _settings.get<double>("preservationRateInit") < 0.000000001){
        std::cout << "Invalid initial settings " << std::endl;
        std::cout << " cannot have <<numberOccurrences>> greater than 0 and " << std::endl;
        std::cout << " <<preservationRateInit>> equal to zero. Check control file " << std::endl;
        exit(0);
    
    }
    
    double updateRatePreservationRate = _settings.get<double>("updateRatePreservationRate");
    
    
    if (_numberOccurrences == 0){
        _hasPaleoData = false;
        _observationTime = _tree->getAge();
        
        if (getTreePtr()->isUltrametric() == false){
            std::cout << "Tree must be ultrametric if no fossil data" << std::endl;
            std::cout << "Exiting...." << std::endl;
            exit(0);
        }
 
        if (updateRatePreservationRate > 0.00000001){
            std::cout << "Attempt to set preservation rate for non-paleo data" << std::endl;
            std::cout << "This parameter will be ignored..." << std::endl;
        }
        
    }else if (_numberOccurrences > 0){
        
        _hasPaleoData = true;
        
        // Observation time of tree:
        _observationTime = _settings.get<double>("observationTime");
        if (_observationTime <= 0){
            _observationTime = _tree->getAge();
        }else if ( _observationTime < _tree->getAge() ){
            std::cout << "WARNING: invalid initial observation time" << std::endl;
            std::cout << "\t... setting observationTime to tree MAX TIME" << std::endl;
            _observationTime = _tree->getAge();
        }
        
    }else if (_numberOccurrences == -1){
        // This is the flag to tell BAMM to expect a formal
        // preservation rate file.
        

        

        _hasPaleoData = true;
        _observationTime = _settings.get<double>("observationTime");
        if (_observationTime <= 0){
            _observationTime = _tree->getAge();
        }else if ( _observationTime < _tree->getAge() ){
            std::cout << "WARNING: invalid initial observation time" << std::endl;
            std::cout << "\t... setting observationTime to tree MAX TIME" << std::endl;
            _observationTime = _tree->getAge();
        }
        
        if (_settings.get("fossilDataFile") == ""){
            std::cout << "Invalid input file for fossil occurrences" << std::endl;
            std::cout << "This must be specified if <<numberOccurrences>> = -1" << std::endl;
            exit(0);
        }else{
        
            getFossilDataFromFile(_settings.get("fossilDataFile"));
 
        }        
        
        
    }else{
        std::cout << "Invalid number of occurrences in controlfile" << std::endl;
        exit(0);
    }
    
    _hasMassExtinctionData = _settings.get<bool>("hasMassExtinctions");
    
    if (_hasMassExtinctionData){
        //std::cout << "file: " << _settings.get("massExtinctionFile");
        getMassExtinctionDataFromFile();
    }
    

}

void SpExModel::getFossilDataFromFile(std::string fileName)
{
    

    std::ifstream inputFile(fileName.c_str());
    
    if (!inputFile) {
        log(Error) << "Could not read data from file "
        << "<<" << fileName << ">>.\n";
        std::exit(1);
    }
    
    log() << "Reading fossil data from file <<" << fileName << ">>.\n";
 
    while (inputFile){
        std::string tempstring;
        getline(inputFile, tempstring, '\t');
        _stagenames.push_back(tempstring);
        getline(inputFile, tempstring, '\t');
        _startTime.push_back(atof(tempstring.c_str()));
        getline(inputFile, tempstring, '\t');
        _endTime.push_back(atof(tempstring.c_str()));
        getline(inputFile, tempstring, '\t');
        _relPresRate.push_back(atof(tempstring.c_str()));
        getline(inputFile, tempstring, '\n');
        _fossilCount.push_back(atof(tempstring.c_str()));

        // this OK?
        if (inputFile.peek() == EOF) {
            break;
        }
    }
    
    inputFile.close();
 
    if (_endTime[(int)_endTime.size()-1] < _observationTime){
        _endTime[(int)_endTime.size()-1] = _observationTime;
        
        std::cout << "WARNING: preservation file should have final time bin" << std::endl;
        std::cout << "consistent with <<observationTime>>" << std::endl;
        std::cout << "Resetting end time of final bin to <<observationTime>> " << std::endl;
 
    }
    
    
}

void SpExModel::getMassExtinctionDataFromFile()
{

    std::string fileName = _settings.get("massExtinctionFile");
 
    std::ifstream inputFile(fileName.c_str());
    
    if (!inputFile) {
        log(Error) << "Could not read data from file "
        << "<<" << fileName << ">>.\n";
        std::exit(1);
    }
    
    log() << "Reading fossil data from file <<" << fileName << ">>.\n";

    while (inputFile){
        std::string tempstring;
        getline(inputFile, tempstring, '\t');
        _massExtinctionTime.push_back(atof(tempstring.c_str()));
        getline(inputFile, tempstring, '\n');
        _massExtinctionIntensity.push_back(atof(tempstring.c_str()));
        
        // this OK?
        if (inputFile.peek() == EOF) {
            break;
        }
    }
    
    inputFile.close();
    
    for (int i = 0; i < (int)_massExtinctionIntensity.size(); i++){
        std::cout << "MassEx Pars" << std::endl;
        std::cout << _massExtinctionIntensity[i] << std::endl;
        std::cout << _massExtinctionTime[i] << std::endl;
    }

    if ((int)_massExtinctionIntensity.size() > 1){
        std::cout << "Only single (terminal) mass extinction currently supported" << std::endl;
        exit(0);
    }

}




double SpExModel::getScaledPreservationRate(double abstime)
{
    double rate = 0.0;
    if (_numberOccurrences > 0){
        rate = _preservationRate;
    }else if (_numberOccurrences < 0){
        for (int i = 0; i < (int)_startTime.size(); i++){
            if (abstime >= _startTime[i] & abstime <= _endTime[i]){
                rate = _relPresRate[i] * _preservationRate;
                break;
            }
        }
 
    }else{
        std::cout << "Invalid _numberOccurrences" << std::endl;
        exit(0);
    }


    return rate;
}

void SpExModel::testPreservationRate(void)
{
    
    for (int i = 0; i < (int)_startTime.size(); i++){
        std::cout << _startTime[i] << "\t" << _endTime[i] << "\t" << _relPresRate[i] << std::endl;
    }
    
    
    double tt = 0.0;
    while (tt < _observationTime){
        std::cout << "Time : " << tt << "\tPresRate: " << getScaledPreservationRate(tt) << std::endl;
        tt += 10;
    }
    
}


// Returns probability that a given lineage alive at start of interval will
// have gone extinct, on account of mass extinctions alone.
double SpExModel::getMassExtinctionPointIntensity(double t_start, double t_end)
{
    double psurvive = 1.0;
    
    for (int i = 0; i < (int)_massExtinctionTime.size(); i++){
        if ((t_start <= _massExtinctionTime[i]) & (t_end >= _massExtinctionTime[i]))
            psurvive *= (1 - _massExtinctionIntensity[i]);
    }
    
    return (1 - psurvive);
    
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


BranchEvent* SpExModel::newBranchEventWithParametersFromSettings(double x)
{
    
    // x is map time
    
    
    double newLam = _settings.get<double>("lambdaInit0");
    double newMu = _settings.get<double>("muInit0");
    double newMuShift = _settings.get<double>("muShift0");
    bool newIsTimeVariable = _prior.generateLambdaIsTimeVariableFromPrior();
    double newLambdaShift = _settings.get<double>("lambdaShift0");
    
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
            
            logLikelihood += (LL + LR);

            // Does not include root node, so it is conditioned
            // on basal speciation event occurring:
            if (node != _tree->getRoot()) {
                logLikelihood  += log(node->getNodeLambda());
                
                // DEBUG
                //std::cout << node << "\tlambda:\t" << node->getNodeLambda();
                //std::cout << "\tEvents:\t" << node->getBranchHistory()->getNumberOfBranchEvents() << std::endl;
                
                node->setDinit(1.0);
            }
        }
    }

    
    if (_hasPaleoData){
        logLikelihood += computePreservationLogProb();  
    }
    
    return logLikelihood;
}


double SpExModel::computeSpExProbBranch(Node* node)
{
    // DEBUG
    //int nodeindex = nodeIndexLookup(node);
    //std::cout << node << "\t" << node->getRandomRightDesc() << "\t" << node->getRandomLeftDesc() << std::endl;
    
    double logLikelihood = 0.0;

    double D0 = node->getDinit();    // Initial speciation probability
    double E0 = node->getEinit();    // Initial extinction probability
    
    

    // 3 scenarios:
    //   i. node is extant tip
    //   ii. node is fossil last occurrence, an unsampled or extinct tip
    //   iii. node is internal. Will now treat separately.
    //  case i and iii can be treated the same
    
    
    // Test if node is extant
    
    // TODO: This test for "extant" vs "non-extant" can be a problem,
    //      depending on numerical error etc. There must be some sort of check.
    //      If lineages that are EXTANT are looped through here,
    //      you will have massively depressed log-likelihoods as it will compute
    //      and add extinction likelihoods for extant lineages.
    //      Problems were observed with simulated trees when the tolerance parameter
    //      was set to 0.00001, as it was flagging many extant taxa as extinct.
    
    bool isExtant = (std::abs(node->getTime() - _observationTime)) < 0.01;

    // TODO: confirm that this is correctly doing the following:
    //        computing the probability of an unobserved lineage
    //         going extinct before the present (or not being sampled)
    //         given a last occurrence time in a tree (eg extinction observed)
    //         and then initializing everything for the downstream calculations
    
    
    if (node->isInternal() == false & isExtant == false){
    // case 1: node is fossil tip
    
        double ddt = _observationTime - node->getTime();
        
        double startTime = node->getBrlen() + ddt;
        double endTime = startTime;
        
        // MASS EXTINCTION INTENSITY
        E0 = getMassExtinctionPointIntensity(node->getTime(), _observationTime);
 
        
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
            
            // TODO: check interval-specific psi values
            // curPsi can be computed once we have interval-specific psi values
            //
            //
            // Old way:
            //double curPsi = _preservationRate;
    
            double atime = node->getAnc()->getTime() + startTime;
            double curPsi = getScaledPreservationRate(atime);
            
            // This (above) definition of atime should work because:
            //    for a given node, startTime = 0 occurs at the ROOTWARD
            //    end of branch starting with the node. Hence, at startTIme = 0
            //    you are at the parent node (node->getAnc()).
            // Then you simply add this time to the absolute age of the parent node
            //    to get the absolute time.
            
            
            double spProb = 0.0;
            double exProb = 0.0;
        
            computeSpExProb(spProb, exProb, curLam, curMu, curPsi, D0, E0, deltaT);

            // exProb cannot exceed max probability of extinction
            //   necessary to avoid overflow/underflow issues
            //   as E0 approaches 1.0
            // Hence we automatically
            //   set any such values such that the computed log-likelihood
            //   is -INFINITY to avoid overflow/underflow at the boundaries
            //   of permitted parameters.
            // Fix added June 13 2015
            // TODO: re-re-recheck this for fossilBAMM
            //  also re-check E0 initializations for
            //  extinct terminal branches
            
            //if (exProb > _extinctionProbMax || spProb > 1.0 ){
            // This line above is incorrect: spProb CAN exceed 1.0
            // if sum of D0 and E0 > 1. I think theoretical bound on
            // value of D(t) is E0 + D0 (check).
            
            if (exProb > _extinctionProbMax ){
                
                return -INFINITY;
            }
            
            E0 = exProb;
            
            endTime = startTime;
            
        }
        
        // Prob that lineage went extinct before present
        // E0 could be the new D0 for the next calculation
        //  however, we will factor this out and start with 1.0.
        
        logLikelihood += log(E0);

        
        // TODO:: the next 2 lines: why were they set? In class Tree,
        //      check that Di and Ei values for tip nodes (extinct or extant) are set correctly.
        //node->setDinit(1.0);
        //node->setEinit(E0); // Why is this getting reset???
        
        D0 = 1.0;
        // current value of E0 can now be passed on for
        // calculation down remainder of branch (eg, the observed segement)
 
    }
    
    // DEBUG
    //_Einitvec[nodeindex] = E0;
    //_Dinitvec[nodeindex] = D0;
    
    
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

        
        // TODO: doublecheck calculations of curPsi for stage-specific rates.
        //double curPsi = _preservationRate;
        
        double atime = node->getAnc()->getTime() + startTime;
        double curPsi = getScaledPreservationRate(atime);
        
        // DEBUG
        //std::cout << "Node: \t" << node << "\tDO\t" << D0 << "\tE0\t" << E0 << std::endl;
        
        
        double spProb = 0.0;
        double exProb = 0.0;

        // Compute speciation and extinction probabilities and store them
        // in spProb and exProb (through reference passing)
        computeSpExProb(spProb, exProb, curLam, curMu, curPsi, D0, E0, deltaT);

        // exProb cannot exceed max probability of extinction
        //   necessary to avoid overflow/underflow issues
        //   as E0 approaches 1.0
        // Hence we automatically
        //   set any such values such that the computed log-likelihood
        //   is -INFINITY to avoid overflow/underflow at the boundaries
        //   of permitted parameters.
        // Fix added June 13 2015
        
        // DEBUG
        //std::cout << "exProb/spProb " << exProb << "\t" << spProb << "\tST/ET\t" << startTime << "\t" << endTime << std::endl;
        
        
        //if (exProb > _extinctionProbMax || spProb > 1.0 ){
        // This line above is incorrect: spProb CAN exceed 1.0
        // if sum of D0 and E0 > 1. I think theoretical bound on
        // value of D(t) is E0 + D0 (check).
        
        if (exProb > _extinctionProbMax ){
        
            return -INFINITY;
        }
        
        logLikelihood += std::log(spProb);
        
        //std::cout << node << "curLam:\t" << curLam << "\tcurMu:\t" << curMu << "\tspProb:\t" << spProb << std::endl;

        //std::cout << node << "E0:\t" << E0 << "\tD0:\t" << D0 << "sprob:\t" << std::log(spProb) << std::endl;
        
        // DEBUG
        //_Dfinalvec[nodeindex] += std::log(spProb);
        
        
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
            // The value E0 at the end of this calculation is passed down through the
            //   tree so one does need to recompute it.
            
            //TODO : doublecheck this section... make sure no variables recycled inappropriately
            
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
                
                // TODO: doublecheck stage-specific preservation rates.
                //double cpsi = _preservationRate;
                
                double atime = node->getAnc()->getTime() + startTime;
                double cpsi = getScaledPreservationRate(atime);
                
                
                double sprob = 0.0;
                double eprob = 0.0;
    
                computeSpExProb(sprob, eprob, clam, cmu, cpsi, (double)1.0, E0, tt);
                
                // DEBUG
                //std::cout << node << "\tEprob:\t" << eprob << std::endl;
                
                E0 = eprob;
                et = st;
                
                // E0 cannot exceed max probability of extinction
                //   necessary to avoid overflow/underflow issues
                //   as E0 approaches 1.0
                if (E0 > _extinctionProbMax){

                    return -INFINITY;
                }
                
            }
            
            
        }

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
    
    // DEBUG
    //_Efinalvec[nodeindex] = E0;
 
    // BUG VERSION
    // Fixed June 13 2015:
    // Did not check for valid E0, exProb, and spProb values
    // during segment calculations on branches - only at ends of branch.
    // This allowed numerical problems causing E0 = 1 and spProb > 1
    // on individual branch segments. This was probably introduced when
    // fossil calculations were introduced. Unlikely to have compromised any
    // analyses as led to easily diagnosable pathologies when E0 = 1 on
    // branch segments
    //
    // BUG VERSION included only this check below - now commented out:
    //if (E0 > _extinctionProbMax) {
    //    return -INFINITY;
    //}
    
    
    
    // *************************************************// 

    // To CONDITION on survival, uncomment the line below;
    // or to NOT condition on survival, comment the line below
 
    // TODO: add conditional here to condition on survival;
    //          but must test for paleo tree...
    
    
    if (parent == _tree->getRoot() & _conditionOnSurvival == true){
 
        logLikelihood -= std::log(1.0 - E0);
    }
    
    // *************************************************//
    
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

    // DEBUG
    //std::cout << "\n\tFF\t" << FF << "\tc1\t" << c1 << "\tc2\t" << c2;
    //std::cout << "\n\tA\t" << A << "\tB\t" << B << "\tX\t" << X << "\tY\t" << Y << "\n" << std::endl;
    
    
    
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
    
    // Prior density on the preservation rate, if paleo data:
    if (_hasPaleoData){
        logPrior += _prior.preservationRatePrior(_preservationRate);
    }

    return logPrior;
}


double SpExModel::computePreservationLogProb()
{
    double logLik = 0.0;

    if (_numberOccurrences > 0){
        
        logLik = (double)_numberOccurrences * std::log(_preservationRate);
        
    }else if (_numberOccurrences < 0){
        
        for (int i = 0; i < (int)_relPresRate.size(); i++){

            if (_fossilCount[i] > 0){
                
                double logprate = std::log(_preservationRate) + std::log(_relPresRate[i]);
                
                logLik += logprate * (double)_fossilCount[i];
            }
            

        }
        
    }else{
        std::cout << "Invalid _numberOccurrences variable in SpExModel::computePreservationLogProb";
        std::cout << std::endl;
        exit(0);
    }
    
    
    return logLik;
}

/************************/
// DEBUG FUNCTIONS BELOW

void SpExModel::printNodeProbs()
{
    int numNodes = _tree->getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        Node* node = postOrderNodes[i];
        std::cout << node << "\t" << node->getName() << "\t";
        std::cout << node->getEinit() << std::endl;
    
    }

}



void SpExModel::initializeDebugVectors()
{
    int numNodes = _tree->getNumberOfNodes();
    for (int i = 0; i < numNodes; i++){
        _Dinitvec.push_back(0.0);
        _Dfinalvec.push_back(0.0);
        _Einitvec.push_back(0.0);
        _Efinalvec.push_back(0.0);
    }

}

int SpExModel::nodeIndexLookup(Node* node)
{
    int goodnode = -1;
    int numNodes = _tree->getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        if (node == postOrderNodes[i]){
            goodnode = i;
        }
    }
    if (goodnode == -1){
        std::cout << "failed to match node." << std::endl;
        exit(0);
    }
    return goodnode;
}

void SpExModel::outputDebugVectors()
{
    
    std::string outname = "debug_LH.txt";
    std::ofstream outStream;
    outStream.open(outname.c_str());
    
    // Write header:
    outStream << "Rdesc,Ldesc,blen,D_init,D_final,E_init,E_final\n";
    //outStream.close();
 
    int numNodes = _tree->getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < numNodes; i++) {
        Node* x = postOrderNodes[i];
        outStream << x->getRandomRightDesc() << ",";
        outStream << x->getRandomLeftDesc() << ",";
        outStream << x->getBrlen() << ",";
        outStream << _Dinitvec[i] << ",";
        outStream << _Dfinalvec[i] << ",";
        outStream << _Einitvec[i] << ",";
        outStream << _Efinalvec[i] << "\n";
        
        // x->getRandomLeftDesc() & x->getRandomRightDesc()
        
    }
    
    outStream.close();
    
    
}





