#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cctype>
#include <algorithm>
#include <cstdlib>

#include "Settings.h"
#include "Log.h"
#include "MatchPathSeparator.h"


Settings::Settings(const std::string& controlFilename,
    const std::vector<UserParameter>& commandLineParameters) :
    _commandLineParameters(commandLineParameters)
{
    readControlFile(controlFilename);

    // Get the model type
    std::string modelType;
    std::vector<UserParameter>::const_iterator it;
    for (it = _userParameters.begin(); it < _userParameters.end(); ++it) {
        if (it->first == "modeltype") {
            modelType = it->second;
            break;
        }
    }
 
    initializeGlobalSettings();
 
    // Initialize specific settings for model type
    if (modelType == "speciationextinction") {
        initializeSpeciationExtinctionSettings();
    } else if (modelType == "trait") {
        initializeTraitSettings();
    } else {
        exitWithErrorInvalidModelType();
    }

    // Re-assign parameters based on user values
    initializeSettingsWithUserValues();

    checkAllSettingsAreUserDefined();
    checkAllOutputFilesAreWriteable();
    
    validateSettings();
    
}


void Settings::readControlFile(const std::string& controlFilename)
{
    if (!fileExists(controlFilename)) {
        exitWithErrorNoControlFile();
    }

    std::ifstream controlStream(controlFilename.c_str());

    while (!controlStream.eof()) {
        std::string line;
        getline(controlStream, line, '\n');

        // Strip whitespace
        line.erase(std::remove_if(line.begin(), line.end(),
            (int(*)(int))isspace), line.end());

        // Strip comments
        std::istringstream lineStringStream(line);
        std::string lineWithoutComments;
        getline(lineStringStream, lineWithoutComments, '#');

        // Skip empty lines
        if (lineWithoutComments.size() == 0) {
            continue;
        }

        // Use second getline to split by '=' characters
        std::istringstream equalsStringStream(lineWithoutComments);
        std::vector<std::string> tokens;
        std::string token;
        while (getline(equalsStringStream, token, '=')) {
            tokens.push_back(token);
        }

        // Ensure input line is valid (two sides to an equal sign)
        if (tokens.size() != 2) {
            exitWithErrorInvalidLine(lineWithoutComments);
        }

        // Store parameter and its value
        _userParameters.push_back(UserParameter(tokens[0], tokens[1]));
    }
}


void Settings::initializeGlobalSettings()
{
    // General
    addParameter("modeltype", "speciationextinction");
    addParameter("treefile", "tree.txt");
    addParameter("sampleFromPriorOnly", "0", NotRequired);
    addParameter("runMCMC", "0");
    addParameter("loadEventData", "0", NotRequired);
    addParameter("eventDataInfile", "event_data_in.txt", NotRequired);
    addParameter("initializeModel", "0");
    addParameter("simulatePriorShifts", "0", NotRequired);
    addParameter("numberOfGenerations", "0");
    addParameter("seed", "-1", NotRequired);
    addParameter("validateEventConfiguration", "0", NotRequired);
    addParameter("checkUltrametric", "1", NotRequired);

    // MCMC tuning
    addParameter("updateEventLocationScale", "0.0");
    addParameter("updateEventRateScale", "0.0");
    addParameter("localGlobalMoveRatio", "0.0");

    // Metropolis-coupled MCMC
    addParameter("numberOfChains", "1", NotRequired);
    addParameter("deltaT", "0.1", NotRequired);
    addParameter("swapPeriod", "1000", NotRequired);
    addParameter("chainSwapFileName", "chain_swap.txt", NotRequired);

    // Priors
    addParameter("poissonRatePrior", "0.0", NotRequired);
    addParameter("expectedNumberOfShifts", "0.0", NotRequired);

    // Output
    addParameter("outName", "", NotRequired);
    addParameter("runInfoFilename", "run_info.txt", NotRequired);
    addParameter("mcmcOutfile", "mcmc_out.txt", NotRequired);
    addParameter("eventDataOutfile", "event_data.txt", NotRequired);
    addParameter("priorOutputFileName", "prior_probs.txt", NotRequired);

    addParameter("branchRatesWriteFreq", "0", NotRequired);
    addParameter("mcmcWriteFreq", "0");
    addParameter("eventDataWriteFreq", "0");

    addParameter("printFreq", "0");
    addParameter("overwrite", "0", NotRequired);
    addParameter("writeMeanBranchLengthTrees", "0", NotRequired);

    addParameter("acceptanceResetFreq", "1000", NotRequired);

    // Parameter update rates
    addParameter("updateRateEventNumber", "0.0");
    addParameter("updateRateEventNumberForBranch", "0.0", NotRequired);
    addParameter("updateRateEventPosition", "0.0");
    addParameter("updateRateEventRate", "0.0");
    addParameter("initialNumberEvents", "0");

    // Other (TODO: Need to add documentation for these)
    addParameter("autotune", "0", NotRequired);
    addParameter("outputAcceptanceInfo", "0", NotRequired);
    addParameter("acceptanceInfoFileName", "acceptance_info.txt", NotRequired);

    // TODO: New params May 30 2014, need documented
    addParameter("maxNumberEvents", "5000", NotRequired);
    addParameter("priorSim_IntervalGenerations", "5000", NotRequired);
    addParameter("fastSimulatePrior_Generations", "5000000", NotRequired);
    addParameter("fastSimulatePrior_SampleFreq", "50", NotRequired);
    addParameter("fastSimulatePriorExperimental", "0", NotRequired);
    addParameter("fastSimulatePrior_BurnIn", "0.05", NotRequired);

    // maxNumberEvents = for fastSimulatePriorExperimental, maximum number of events...
    // priorSims_intervalGenerations = number of generations per model pair
    // fastSimulatePrior_Generations = sampling gens for simulation of prior
    // fastSimulatePrior_SampleFreq = sampling frequency for fastSimulation of prior
    // fastSimulatePrior_BurnIn = fraction of samples to be discarded as burnin.

    /********************************************************/
    // Parameters for fossilBAMM
    // TODO: this should NOT be a global set of parameters.
    // But it will not run trait BAMM unless they are set here.
    // Problem is with
    
    addParameter("preservationRateInit", "0", NotRequired);
    addParameter("observationTime", "-1", NotRequired);
    addParameter("numberOccurrences", "0", NotRequired);
    addParameter("updateRatePreservationRate", "-1", NotRequired);
    addParameter("updatePreservationRateScale", "1.0", NotRequired);
    addParameter("preservationRatePrior", "1.0", NotRequired);
    
    /*********************************************************/
    
    
}


void Settings::initializeSpeciationExtinctionSettings()
{
    // General
    addParameter("useGlobalSamplingProbability", "1");
    addParameter("globalSamplingFraction", "0.0");
    addParameter("sampleProbsFilename", "sample_probs.txt", NotRequired);

    // MCMC tuning
    addParameter("updateLambdaInitScale", "0.0");
    addParameter("updateMuInitScale", "0.0");
    addParameter("updateLambdaShiftScale", "0.0");
    addParameter("updateMuShiftScale", "0.0", NotRequired);
    addParameter("minCladeSizeForShift", "1", NotRequired);

    // Starting parameters
    addParameter("lambdaInit0", "0.0");
    addParameter("lambdaShift0", "0.0");
    addParameter("muInit0", "0.0");
    addParameter("muShift0", "0.0", NotRequired);

    // Priors
    addParameter("lambdaInitPrior", "0.0");
    addParameter("lambdaShiftPrior", "0.0");
    addParameter("muInitPrior", "0.0");
    addParameter("muShiftPrior", "1.0", NotRequired);
    addParameter("lambdaInitRootPrior", "-1.0", NotRequired);
    addParameter("lambdaShiftRootPrior", "-1.0", NotRequired);
    addParameter("muInitRootPrior", "-1.0", NotRequired);
    addParameter("muShiftRootPrior", "-1.0", NotRequired);
    addParameter("lambdaIsTimeVariablePrior", "0.5");
	addParameter("segLength", "0.0");

    // Output
    addParameter("lambdaOutfile", "lambda_rates.txt", NotRequired, Deprecated);
    addParameter("muOutfile", "mu_rates.txt", NotRequired, Deprecated);

    // Parameter update rates
    addParameter("updateRateLambda0", "0.0");
    addParameter("updateRateLambdaShift", "0.0");
    addParameter("updateRateMu0", "0.0");
    addParameter("updateRateMuShift", "0.0", NotRequired);
    addParameter("updateRateLambdaTimeMode", "0.0");

    // Maximum value of extinction probability on branch that will be tolerated:
    // to avoid numerical overflow issues (especially rounding to 1)
    addParameter("extinctionProbMax", "0.9999", NotRequired);
    
    addParameter("conditionOnSurvival", "-1", NotRequired);
    addParameter("alwaysRecomputeE0", "0", NotRequired);
    
    addParameter("combineExtinctionAtNodes", "if_different", NotRequired);
    
    
    /********************************************************/
    // Parameters for fossilBAMM
    
    /*
    addParameter("preservationRateInit", "0", NotRequired);
    addParameter("observationTime", "-1", NotRequired);
    addParameter("numberOccurrences", "0", NotRequired);
    addParameter("updateRatePreservationRate", "-1", NotRequired);
    addParameter("updatePreservationRateScale", "1.0", NotRequired);
    addParameter("preservationRatePrior", "1.0", NotRequired);
    */
     
    /*********************************************************/
    
}


void Settings::initializeTraitSettings()
{
    // General
    addParameter("traitfile", "traits.txt");

    // MCMC tuning
    addParameter("updateBetaInitScale", "0.0");
    addParameter("updateNodeStateScale", "0.0");
    addParameter("updateBetaShiftScale", "0.0");

    // Starting parameters
    addParameter("betaInit", "0.0");
    addParameter("betaShiftInit", "0.0");

    // Priors
    addParameter("betaInitPrior", "0.0");
    addParameter("betaShiftPrior", "0.0");
    addParameter("betaInitRootPrior", "-1.0", NotRequired);
    addParameter("betaShiftRootPrior", "-1.0", NotRequired);
    addParameter("betaIsTimeVariablePrior", "0.5");
    addParameter("useObservedMinMaxAsTraitPriors", "1");
    addParameter("traitPriorMin", "0.0", NotRequired);
    addParameter("traitPriorMax", "0.0", NotRequired);

    // Output
    addParameter("betaOutfile", "beta_rates.txt", NotRequired, Deprecated);
    addParameter("nodeStateOutfile", "node_state.txt", NotRequired);
    addParameter("nodeStateWriteFreq", "0", NotRequired);

    // Parameter update rates
    addParameter("updateRateBeta0", "0.0");
    addParameter("updateRateBetaShift", "0.0");
    addParameter("updateRateNodeState", "0.0");
    addParameter("updateRateBetaTimeMode", "0.0");
}


void Settings::addParameter(const std::string& name, const std::string& value,
    UserDefinedStatus userDefined, DeprecationStatus deprecated)
{
    _parameters.insert(Parameter(name,
        SettingsParameter(name, value, userDefined, deprecated)));
}



void Settings::initializeSettingsWithUserValues()
{
    std::vector<std::string> paramsNotFound;

    std::vector<UserParameter>::const_iterator userParamIt;
    for (userParamIt = _userParameters.begin();
        userParamIt != _userParameters.end(); ++userParamIt) {

        // Find the matching parameter to the user-specified parameter
        ParameterMap::iterator paramIt = _parameters.find((userParamIt->first));

        // If found, set the value of the parameter to the user's
        if (paramIt != _parameters.end()) {
            // Parameter should not be deprecated
            if ((paramIt->second).isDeprecated()) {
                exitWithErrorParameterIsDeprecated(paramIt->first);
            }

            // Parameter should not already be user-defined
            if ((paramIt->second).isUserDefined()) {
                exitWithErrorDuplicateParameter(paramIt->first);
            }

            (paramIt->second).setStringValue(userParamIt->second);
        } else {
            paramsNotFound.push_back(userParamIt->first);
        }
    }

    // Handle command-line arguments
    std::vector<UserParameter>::const_iterator cmdLineParamIt;
    for (cmdLineParamIt = _commandLineParameters.begin();
        cmdLineParamIt != _commandLineParameters.end(); ++cmdLineParamIt) {
        // Find the command-line parameter in known parameter list
        ParameterMap::iterator paramIt =
            _parameters.find((cmdLineParamIt->first));

        // If found, replace current parameter value with command-line value
        if (paramIt != _parameters.end()) {
            (paramIt->second).setStringValue(cmdLineParamIt->second);
        } else {
            log(Error) << "Command-line parameter " << cmdLineParamIt->first
                << " is not a known parameter.\n";
            std::exit(1);
        }
    }

    attachPrefixToOutputFiles();

    if (paramsNotFound.size() > 0) {
        exitWithErrorParametersNotFound(paramsNotFound);
    }
}


void Settings::attachPrefixToOutputFiles()
{
    // Get the prefix string value (in outName parameter)
    ParameterMap::const_iterator it = _parameters.find("outName");
    std::string prefix;
    if (it != _parameters.end()) {
        prefix = (it->second).value<std::string>();
    }

    // Create an array of the parameters that need to be prefixed
    std::string paramsToPrefix[NumberOfParamsToPrefix] =
        { "runInfoFilename",
          "mcmcOutfile",
          "eventDataOutfile",
          "nodeStateOutfile",
          "priorOutputFileName",
          "acceptanceInfoFileName",
          "chainSwapFileName",
          "lambdaOutfile",
          "muOutfile",
          "betaOutfile" };

    // Attach the prefix to each parameter
    ParameterMap::iterator paramIt;
    for (size_t i = 0; i < NumberOfParamsToPrefix; i++) {
        paramIt = _parameters.find(paramsToPrefix[i]);
        if (paramIt != _parameters.end()) {
            const std::string& param = (paramIt->second).value<std::string>();
            (paramIt->second).setStringValue(attachPrefix(prefix, param));
        }
    }
}


std::string Settings::attachPrefix
  (const std::string& prefix, const std::string& path) const
{
    if (prefix == "") {
        return path;
    }

    const std::string& dir = extractDir(path);
    const std::string& fileName = extractFileName(path);

    return dir + prefix + "_" + fileName;
}


std::string Settings::extractDir(const std::string& path) const
{
    // Starts from the beginning, stops when it finds the last path separator
    return std::string(path.begin(), std::find_if(path.rbegin(), path.rend(),
        MatchPathSeparator()).base());
}


std::string Settings::extractFileName(const std::string& path) const
{
    // Starts from the end and stops when it finds a path separator
    return std::string(std::find_if(path.rbegin(), path.rend(),
        MatchPathSeparator()).base(), path.end());
}


void Settings::checkAllSettingsAreUserDefined() const
{
    ParameterMap::const_iterator it;
    for (it = _parameters.begin(); it != _parameters.end(); ++it) {
        const SettingsParameter& parameter = it->second;
        if (!parameter.isDeprecated()) {
            if (!parameter.isUserDefined() && parameter.mustBeUserDefined()) {
                exitWithErrorUndefinedParameter(parameter.name());
            }
        }
    }
}


void Settings::checkAllOutputFilesAreWriteable() const
{
    if (!get<bool>("overwrite")) {
        if (anyOutputFileExists()) {
            exitWithErrorOutputFileExists();
        }
    }
}


bool Settings::anyOutputFileExists() const
{
    // Global output files
    if (fileExists(get("runInfoFilename")) ||
        fileExists(get("mcmcOutfile"))     ||
        fileExists(get("eventDataOutfile"))) {
        return true;
    }

    if (get<bool>("writeMeanBranchLengthTrees")) {
        // Speciation/extinction output files
        if (get("modeltype") == "speciationextinction") {
            if (fileExists(get("lambdaOutfile")) ||
                fileExists(get("muOutfile"))) {
                return true;
            }

        // Trait output files
        } else if (get("modeltype") == "trait") {
            if (fileExists(get("betaOutfile"))) {
                return true;
            }
        }
    }

    return false;
}



bool Settings::fileExists(const std::string& filename) const
{
    std::ifstream inFile(filename.c_str());
    return inFile.good();
}


void Settings::printCurrentSettings(std::ostream& out) const
{
    int ppw = 29;

    log(Message, out) << "Current parameter settings:\n";

    ParameterMap::const_iterator it;
    for (it = _parameters.begin(); it != _parameters.end(); ++it) {
        log(Message, out) << std::right << std::setw(ppw) <<
            it->first << "    " << (it->second).value<std::string>() << "\n";
    }

    log(Message, out) << std::flush;
}


void Settings::exitWithErrorNoControlFile() const
{
    log(Error) << "Specified control file does not exist.\n"
               << "Check that the file is in the specified location.\n";
    std::exit(1);
}


void Settings::exitWithErrorInvalidLine(const std::string& line) const
{
    log(Error) << "Invalid input line in control file.\n"
               << "Problematic line includes <<" << line << ">>\n";
    std::exit(1);
}


void Settings::exitWithErrorUndefinedParameter(const std::string& name) const
{
    log(Error) << "Parameter " << name << " is undefined.\n"
               << "Fix by giving the parameter a value in the control file.\n";
    std::exit(1);
}


void Settings::exitWithErrorInvalidModelType() const
{
    log(Error) << "Invalid type of analysis.\n"
               << "Fix by setting modeltype as speciationextinction or trait\n";
    std::exit(1);
}


void Settings::exitWithErrorParametersNotFound
    (const std::vector<std::string>& paramsNotFound) const
{
    log(Error) << "One or more parameters from the control file does not\n"
        << "correspond to valid model parameters. Make sure that you are\n"
        << "running the correct version of BAMM or check that the following\n"
        << "parameters are spelled correctly:\n\n";

    std::vector<std::string>::const_iterator it;
    for (it = paramsNotFound.begin(); it != paramsNotFound.end(); ++it) {
        log() << std::setw(30) << *it << std::endl;
    }

    std::exit(1);
}

void Settings::exitWithErrorParameterIsDeprecated
    (const std::string& param) const
{
    log(Error) << "Parameter " << param << " has been deprecated. Fix by\n"
        << "removing this parameter or use the appropriate version of BAMM.\n";
    std::exit(1);
}


void Settings::exitWithErrorDuplicateParameter(const std::string& param) const
{
    log(Error) << "Duplicate parameter " << param << ".\n"
               << "Fix by removing duplicate parameter in control file.\n";
    std::exit(1);
}


void Settings::exitWithErrorOutputFileExists() const
{
    log(Error) << "Analysis is set to not overwrite files.\n"
               << "Fix by removing or renaming output file(s),\n"
               << "or set \"overwrite = 1\" in the control file.\n";
    std::exit(1);
}



void Settings::set(const std::string& name, const std::string& value)
{
    ParameterMap::iterator it = _parameters.find(name);
    if (it != _parameters.end()) {
        
        (it->second).setStringValue(value);
        
    }
}

void Settings::validateSettings(void)
{
    //std::cout << "Initial P/E\t" << this->get("poissonRatePrior");
    //std::cout << "\t" << this->get("expectedNumberOfShifts") << std::endl;
    
    double minval = 0.000001;
    double poisson = this->get<double>("poissonRatePrior");
    double expected = this->get<double>("expectedNumberOfShifts");
    
    if (poisson < minval & expected < minval){
        std::cout << "You must specify either:\n";
        std::cout << "\tpoissonRatePrior = <value>, or\n";
        std::cout << "\texpectedNumberOfShifts = <value>\n";
        exit(0);
    }
    
    if (poisson < minval & expected >= minval){
        poisson = 1 / expected;
        
        std::ostringstream s;
        s << poisson;
        this->set("poissonRatePrior", s.str());
    }
 
    //std::cout << "End P/E\t" << this->get("poissonRatePrior");
    //std::cout << "\t" << this->get("expectedNumberOfShifts") << std::endl;
    

}









