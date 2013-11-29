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


Settings::Settings(const std::string& controlFilename)
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
    addParameter("modeltype",                    "speciationextinction");
    addParameter("treefile",                     "tree.txt");
    addParameter("sampleFromPriorOnly",          "0");
    addParameter("runMCMC",                      "0");
    addParameter("loadEventData",                "0");
    addParameter("eventDataInfile",              "event_data_in.txt");
    addParameter("initializeModel",              "0");
    addParameter("numberGenerations",            "0");
    addParameter("seed",                         "-1", false);

    // MCMC tuning
    addParameter("updateEventLocationScale",     "0.0");
    addParameter("updateEventRateScale",         "0.0");
    addParameter("localGlobalMoveRatio",         "0.0");

    // Priors
    addParameter("poissonRatePrior",             "0.0");

    // Output
    addParameter("outName",                      "", false);
    addParameter("runInfoFilename",              "run_info.txt", false);
    addParameter("mcmcOutfile",                  "mcmc_out.txt", false);
    addParameter("eventDataOutfile",             "event_data.txt", false);

    addParameter("branchRatesWriteFreq",         "0", false);
    addParameter("mcmcWriteFreq",                "0");
    addParameter("eventDataWriteFreq",           "0");
    
    addParameter("acceptWriteFreq",              "0");
    addParameter("printFreq",                    "0");
    addParameter("overwrite",                    "0", false);
    addParameter("writeMeanBranchLengthTrees",   "0", false);

    // Parameter update rates
    addParameter("updateRateEventNumber",        "0.0");
    addParameter("updateRateEventPosition",      "0.0");
    addParameter("updateRateEventRate",          "0.0");
    addParameter("initialNumberEvents",          "0");
 
    // Other (TODO: Need to add documentation for these)
    addParameter("autotune",                     "0", false);
}


void Settings::initializeSpeciationExtinctionSettings()
{
    // General
    addParameter("useGlobalSamplingProbability", "1");
    addParameter("globalSamplingFraction",       "0.0");
    addParameter("sampleProbsFilename",          "sample_probs.txt", false);

    // MCMC tuning
    addParameter("updateLambdaInitScale",        "0.0");
    addParameter("updateMuInitScale",            "0.0");
    addParameter("updateLambdaShiftScale",       "0.0");
    addParameter("updateMuShiftScale",           "0.0", false);
    addParameter("minCladeSizeForShift",         "1", false);

    // Starting parameters
    addParameter("lambdaInit0",                  "0.0");
    addParameter("lambdaShift0",                 "0.0");
    addParameter("muInit0",                      "0.0");
    addParameter("muShift0",                     "0.0", false);

    // Priors
    addParameter("lambdaInitPrior",              "0.0");
    addParameter("lambdaShiftPrior",             "0.0");
    addParameter("muInitPrior",                  "0.0");
    addParameter("muShiftPrior",                 "1.0", false);
	addParameter("segLength",                    "0.0");

    // Output
    addParameter("lambdaOutfile",                "lambda_rates.txt", false);
    addParameter("muOutfile",                    "mu_rates.txt", false);

    // Parameter update rates
    addParameter("updateRateLambda0",            "0.0");
    addParameter("updateRateLambdaShift",        "0.0");
    addParameter("updateRateMu0",                "0.0");
    addParameter("updateRateMuShift",            "0.0", false);
}


void Settings::initializeTraitSettings()
{
    // General
    addParameter("traitfile",                      "traits.txt");

    // MCMC tuning
    addParameter("updateBetaScale",                "0.0");
    addParameter("updateNodeStateScale",           "0.0");
    addParameter("updateBetaShiftScale",           "0.0");

    // Starting parameters
    addParameter("betaInit",                       "0.0");
    addParameter("betaShiftInit",                  "0.0");

    // Priors
    addParameter("betaInitPrior",                  "0.0");
    addParameter("betaShiftPrior",                 "0.0");
    addParameter("useObservedMinMaxAsTraitPriors", "1");
    addParameter("traitPriorMin",                  "0.0");
    addParameter("traitPriorMax","0.0");

    // Output
    addParameter("betaOutfile",                    "beta_rates.txt", false);

    // Parameter update rates
    addParameter("updateRateBeta0",                "0.0");
    addParameter("updateRateBetaShift",            "0.0");
    addParameter("updateRateNodeState",            "0.0");
}


void Settings::addParameter(const std::string& name, const std::string& value,
    bool mustBeUserDefined)
{
    _parameters.insert(Parameter(name,
        SettingsParameter(name, value, mustBeUserDefined)));
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
            // Parameter should not already be user-defined
            if ((paramIt->second).isUserDefined()) {
                exitWithErrorDuplicateParameter(paramIt->first);
            } else {
                (paramIt->second).setStringValue(userParamIt->second);
            }
        } else {
            paramsNotFound.push_back(userParamIt->first);
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
  (const std::string& prefix, const std::string& str) const
{
    return (prefix != "") ? (prefix + "_" + str) : (str);
}


void Settings::checkAllSettingsAreUserDefined() const
{
    ParameterMap::const_iterator it;
    for (it = _parameters.begin(); it != _parameters.end(); ++it) {
        const SettingsParameter& parameter = it->second;
        if (!parameter.isUserDefined() && parameter.mustBeUserDefined()) {
            exitWithErrorUndefinedParameter(parameter.name());
        }
    }
}


void Settings::checkAllOutputFilesAreWriteable() const
{
    if (!getOverwrite()) {
        if (anyOutputFileExists()) {
            exitWithErrorOutputFileExists();
        }
    }
}


bool Settings::anyOutputFileExists() const
{
    // Global output files
    if (fileExists(getRunInfoFilename()) ||
        fileExists(getMCMCoutfile())     ||
        fileExists(getEventDataOutfile())) {
        return true;
    }

    if (getWriteMeanBranchLengthTrees()) {
        // Speciation/extinction output files
        if (getModeltype() == "speciationextinction") {
            if (fileExists(getLambdaOutfile()) || fileExists(getMuOutfile())) {
                return true;
            }

        // Trait output files
        } else if (getModeltype() == "trait") {
            if (fileExists(getBetaOutfile())) {
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
               << "Fix by giving the parameter a value in the control file\n";
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
    log(Error) << "One or more parameters from control file\n"
               << "does not correspond to valid model parameters.\n"
               << "Fix by checking the following to see if they are\n"
               << "specified (or spelled) correctly:\n\n";

    std::vector<std::string>::const_iterator it;
    for (it = paramsNotFound.begin(); it != paramsNotFound.end(); ++it) {
        log() << std::setw(30) << *it << std::endl;
    }

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
