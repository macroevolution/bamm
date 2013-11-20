#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cctype>
#include <algorithm>
#include <cstdlib>

#include "Settings.h"


void Settings::initializeSettingsDevel(std::string controlFilename)
{
    std::ifstream infile(controlFilename.c_str());
    if (!infile.good()) {
        std::cout << "Control filename invalid" << std::endl;
        std::cout << "Exiting." << std::endl;
        throw;
    }
    
    std::cout << "Reading control file <<" << controlFilename.c_str() << ">>" << std::endl;
    
    std::string s1, s2, s_nocomment;
    
    while (infile) {
        getline(infile, s1, '\n');
        
        // strip whitespace out of tempstring:
        //      both spaces and tabs:
        
        // What is the int(*)(int) doing????
        s1.erase(std::remove_if(s1.begin(), s1.end(), (int(*)(int))isspace), s1.end());
        
        std::istringstream sx(s1);       
        getline(sx, s_nocomment, '#');
        
        // Only add if has size > 0 (gets rid of empty lines)
        if (s_nocomment.size() > 0) {
            std::vector<std::string> tmpstr;
            
            // Now use second getline to split by '=' characters:

            std::istringstream stemp(s_nocomment);

            while (getline(stemp, s2, '=')) {
                tmpstr.push_back(s2);        
            }

            if (tmpstr.size() == 2) {
                _varName.push_back(tmpstr[0]);
                _varValue.push_back(tmpstr[1]);
            } else {
                std::cout << "Invalid size of input line in control file" << std::endl;
                std::cout << " Problematic line includes <<" << s_nocomment << ">>" << std::endl;
                std::cout << "Terminating run\n" << std::endl;
                std::exit(0);            
            }
        }
        
        if (infile.peek() == EOF) {
            break;
        }
    }

    std::string modelType;
    
    for (std::vector<std::string>::size_type i = 0; i < _varName.size(); i++) {
        if (_varName[i] == "modeltype") {
            modelType = _varValue[i];
        }
    }

    initializeGlobalSettings();
    
    if (modelType == "speciationextinction") {
        initializeSpeciationExtinctionSettings();
    } else if (modelType == "trait") {
        initializeTraitSettings();
    } else {
        std::cout << "Invalid type of analysis" << std::endl;
        std::cout << "Options: speciationextinction or trait" << std::endl;
        std::exit(1);
    }

    initializeSettingsWithUserValues();
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
    addParameter("runInfoFilename",              "run_info.txt");
    addParameter("mcmcOutfile",                  "mcmc_out.txt");
    addParameter("eventDataOutfile",             "event_data.txt");

    // TODO: deprecate acceptrateOutfile
    addParameter("acceptrateOutfile",            "mcmc_accept.txt", false);
    addParameter("branchRatesWriteFreq",         "0");
    addParameter("mcmcWriteFreq",                "0");
    addParameter("eventDataWriteFreq",           "0");
    
    addParameter("acceptWriteFreq",              "0");
    addParameter("printFreq",                    "0");
    addParameter("overwrite",                    "0", false);

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
    addParameter("sampleProbsFilename",          "sample_probs.txt");

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
    addParameter("lambdaOutfile",                "lambda_rates.txt");
    addParameter("muOutfile",                    "mu_rates.txt");
    addParameter("lambdaNodeOutfile",            "node_lambda.txt", false);

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
    addParameter("betaOutfile",                    "beta_rates.txt");
    addParameter("nodeStateOutfile",               "node_states.txt", false);

    // Parameter update rates
    addParameter("updateRateBeta0",                "0.0");
    addParameter("updateRateBetaShift",            "0.0");
    addParameter("updateRateNodeState",            "0.0");
}


void Settings::addParameter(const std::string& name, const std::string& value,
    bool mustBeUserDefined)
{
    _parameters.insert(std::pair<std::string, SettingsParameter>
        (name, SettingsParameter(name, value, mustBeUserDefined)));
}


void Settings::initializeSettingsWithUserValues()
{
    std::vector<std::string> paramsNotFound;

    std::vector<std::string>::size_type i;
    for (i = 0; i < _varName.size(); i++) {
        ParameterMap::iterator it = _parameters.find(_varName[i]);
        if (it != _parameters.end()) {
            assertNotUserDefined(it->second);
            (it->second).setStringValue(_varValue[i]);
        } else {
            // Parameter not found:
            //      add to list of potentially bad/misspelled params
            //      and print for user.
            paramsNotFound.push_back(_varName[i]);
        }
    }

    std::cout << "Read a total of <<" << _varName.size() <<
         ">> parameter settings from control file\n";

    if (paramsNotFound.size() > 0) {
        std::cout << std::endl << "********************************\n";
        std::cout << "ERROR: one or more parameters from control file\n";
        std::cout << "do not correspond to valid model parameters.\n";
        std::cout << "Fix by checking the following to see if they are\n";
        std::cout << "specified (or spelled) correctly:\n\n";

        std::vector<std::string>::const_iterator it;
        for (it = paramsNotFound.begin(); it != paramsNotFound.end(); ++it) {
            std::cout << std::setw(30) << *it << std::endl;
        }

        std::cout << "\n********************************\n\n";
        std::cout << "Execution of BAMM terminated...\n";
        std::exit(1);
    }
}


void Settings::assertNotUserDefined(const SettingsParameter& parameter) const
{
    if (parameter.isUserDefined()) {
        std::cout << "ERROR: Duplicate parameter " << parameter.name() << ".\n";
        std::cout << "Fix by removing duplicate parameter in control file.\n";
        std::exit(1);
    }
}


void Settings::checkSettingsAreUserDefined() const
{
    ParameterMap::const_iterator it;
    for (it = _parameters.begin(); it != _parameters.end(); ++it) {
        const SettingsParameter& parameter = it->second;
        if (!parameter.isUserDefined() && parameter.mustBeUserDefined()) {
            exitWithErrorUndefinedParameter(parameter.name());
        }
    }
}


std::string Settings::attachPrefix(const std::string& str) const
{
    ParameterMap::const_iterator it = _parameters.find("outName");

    std::string prefix;
    if (it != _parameters.end()) {
        prefix = (it->second).value<std::string>();
    }

    if (prefix == "") {
        return str;
    } else {
        return prefix + "_" + str;
    }
}


void Settings::exitWithErrorUndefinedParameter(const std::string& name) const
{
    std::cout << "ERROR: Parameter " << name << " is undefined.\n";
    std::cout << "Fix by assigning the parameter a value in the control file\n";
    std::exit(1);
}


void Settings::printCurrentSettings(std::ostream& out) const
{
    int ppw = 29;

    out << "*****************************************************\n";
    out << "Current parameter settings:\n";

    ParameterMap::const_iterator it;
    for (it = _parameters.begin(); it != _parameters.end(); ++it) {
        out << std::right << std::setw(ppw) <<
            it->first << "\t\t" << (it->second).value<std::string>() << "\n";
    }

    out << "\n*****************************************************\n";
}


void Settings::parseCommandLineInput(int argc, std::vector<std::string>& instrings)
{
    std::vector<std::string> badFlags;

    for (std::vector<std::string>::size_type i = 0; i < (std::vector<std::string>::size_type)argc; i++) {
        if (instrings[i] == "-control") {
            if (instrings.size() == ( i - 1)) {
                std::cout << "Error: no controfile specified" << std::endl;
                std::cout << "Exiting\n\n" << std::endl;
            }
            std::string controlfile = instrings[i + 1];
            std::cout << "\n\nInitializing BAMM using control file <<" << controlfile << ">>" <<
                 std::endl;
            std::ifstream instream(controlfile.c_str());
            if (!instream) {
                std::cout << "File not found error: cannot locate control file in\n";
                std::cout <<  "specified directory. Exiting BAMM" << std::endl << std::endl;
                std::exit(1);
            } else {
                initializeSettingsDevel(controlfile);
            }
        } else {
            if (i != 0) {
                // If i == 0, we just assume that this is either (i) the path to the executable
                //  or (ii), the commands used to invoke it.
                //  May depend on system/compiler, but first argument should
                //  not be a valid flag...
                badFlags.push_back(instrings[i]);
            }
        }
    }
}
