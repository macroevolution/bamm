#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "SettingsParameter.h"

typedef std::pair<std::string, SettingsParameter> Parameter;
typedef std::map<std::string, SettingsParameter> ParameterMap;

typedef std::pair<std::string, std::string> UserParameter;


class Settings
{

private:

    void readControlFile(const std::string& controlFilename);

    void initializeGlobalSettings();
    void initializeSpeciationExtinctionSettings();
    void initializeTraitSettings();
    void initializeSettingsWithUserValues();

    void checkAllSettingsAreUserDefined() const;
    void checkAllOutputFilesAreWriteable() const;

    void assertNotUserDefined(const SettingsParameter& parameter) const;
    void addParameter(const std::string& name, const std::string& value,
        UserDefinedStatus userDefined = Required,
        DeprecationStatus deprecated = NotDeprecated);

    void attachPrefixToOutputFiles();
    std::string attachPrefix
        (const std::string& prefix, const std::string& str) const;
    std::string extractDir(const std::string& path) const;
    std::string extractFileName(const std::string& path) const;

    bool anyOutputFileExists() const;
    bool fileExists(const std::string& filename) const;

    void exitWithErrorNoControlFile() const;
    void exitWithErrorInvalidLine(const std::string& line) const;
    void exitWithErrorUndefinedParameter(const std::string& name) const;
    void exitWithErrorInvalidModelType() const;
    void exitWithErrorParametersNotFound
        (const std::vector<std::string>& paramsNotFound) const;
    void exitWithErrorParameterIsDeprecated(const std::string& param) const;
    void exitWithErrorDuplicateParameter(const std::string& param) const;
    void exitWithErrorOutputFileExists() const;

    static const size_t NumberOfParamsToPrefix = 8;

    // Parameters that settings knows about
    ParameterMap _parameters;
    
    // Parameters read from the control file
    std::vector<UserParameter> _userParameters;

    // Parameters read from the command line
    std::vector<UserParameter> _commandLineParameters;

public:

    Settings(const std::string& controlFilename,
        const std::vector<UserParameter>& commandLineParameters);

    void printCurrentSettings(std::ostream& out = std::cout) const;

    std::string getRunInfoFilename() const;

    std::string getModeltype() const;
	
    // Functions to access parameters for MCMC/Model/Other
    // from Class Settings object:
    bool getUseGlobalSamplingProbability() const;
    bool getSampleFromPriorOnly() const;
    bool getRunMCMC() const;
    bool getInitializeModel() const;
	bool getAutotune() const;

    bool getOutputAcceptanceInfo() const;
    std::string getAcceptanceInfoFileName() const;
    
    bool getSimulatePriorShifts() const;

    // Load previous settings?
    bool getLoadEventData() const;
    std::string getEventDataInfile() const;

    // Sampling probabilities:
    std::string getTreeFilename() const;
    std::string getSampleProbsFilename() const;
    double getGlobalSamplingFraction() const;

    // Class Model parameters:
    double getUpdateLambdaInitScale() const;
    double getUpdateMuInitScale() const;
    double getUpdateLambdaShiftScale() const;
    double getUpdateMuShiftScale() const;
    double getLambdaInit0() const;
    double getLambdaShift0() const;
    double getMuInit0() const;
	double getMuShift0() const;
	double getUpdateEventLocationScale() const;

	//    double getMeanSpeciationLengthFraction();
    double getUpdateEventRateScale() const;
    double getLocalGlobalMoveRatio() const;
    double getPoissonRatePrior() const;
    double getLambdaInitPrior() const;
    double getLambdaShiftPrior() const;
    double getMuInitPrior() const;
    double getMuShiftPrior() const;
    double getLambdaInitRootPrior() const;
    double getLambdaShiftRootPrior() const;
    double getMuInitRootPrior() const;
    double getMuShiftRootPrior() const;
    double getSegLength() const;

    int getMinCladeSizeForShift() const;
    long int getSeed() const;
    bool getOverwrite() const;
    bool getWriteMeanBranchLengthTrees() const;
	
	// Functions to set MCMC operators for autotuning:
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);
	void setUpdateEventRateScale(double x);
	void setUpdateEventLocationScale(double x);

    // Class MCMC parameters:
    std::string getMCMCoutfile() const;
    std::string getEventDataOutfile() const;
    std::string getPriorOutputFileName() const;
    std::string getLambdaOutfile() const;
    std::string getMuOutfile() const;

    int getBranchRatesWriteFreq() const;
    int getEventDataWriteFreq() const;
    int getMCMCwriteFreq() const;
    int getAcceptWriteFreq() const;
    int getPrintFreq() const;
    int getNGENS() const;

    // Class MCMC parameter update weights:
    double getUpdateRateEventNumber() const;
    double getUpdateRateEventPosition() const;
    double getUpdateRateEventRate() const;
    double getUpdateRateLambda0() const;
    double getUpdateRateLambdaShift() const;
    double getUpdateRateMu0() const;
    double getUpdateRateMuShift() const;

    // Other:
    int getInitialNumberEvents() const;
    double getExtinctionProbMax() const;

    /* Parameters specific to trait evolution module */

    std::string getTraitFile() const;

    double getUpdateBetaInitScale() const;
    double getUpdateNodeStateScale() const;
    double getUpdateBetaShiftScale() const;
    double getBetaInit() const;
    double getBetaShiftInit() const;
    double getBetaInitPrior() const;
    double getBetaShiftPrior() const;
    double getBetaInitRootPrior() const;
    double getBetaShiftRootPrior() const;
    double getTraitPriorMin() const;
    double getTraitPriorMax() const;
    double getUpdateRateBeta0() const;
    double getUpdateRateBetaShift() const;
    double getUpdateRateNodeState() const;
    bool   getUseObservedMinMaxAsTraitPriors() const;
    void   setTraitPriorMin(double x);
    void   setTraitPriorMax(double x);
    std::string getBetaOutfile() const;
};


inline std::string Settings::getRunInfoFilename() const
{
    return _parameters.at("runInfoFilename").value<std::string>();
}


inline bool Settings::getAutotune() const
{
	return _parameters.at("autotune").value<bool>();
}


inline bool Settings::getOutputAcceptanceInfo() const
{
    return _parameters.at("outputAcceptanceInfo").value<bool>();
}


inline std::string Settings::getAcceptanceInfoFileName() const
{
    return _parameters.at("acceptanceInfoFileName").value<std::string>();
}


inline bool Settings::getUseGlobalSamplingProbability() const
{
    return _parameters.at("useGlobalSamplingProbability").value<bool>();
}


inline bool Settings::getSampleFromPriorOnly() const
{
    return _parameters.at("sampleFromPriorOnly").value<bool>();
}

inline bool Settings::getSimulatePriorShifts() const
{
    return _parameters.at("simulatePriorShifts").value<bool>();
}


inline bool Settings::getRunMCMC() const
{
    return _parameters.at("runMCMC").value<bool>();
}


inline bool Settings::getInitializeModel() const
{
    return _parameters.at("initializeModel").value<bool>();
}


inline bool Settings::getLoadEventData() const
{
    return _parameters.at("loadEventData").value<bool>();
}


inline std::string Settings::getEventDataInfile() const
{
    return _parameters.at("eventDataInfile").value<std::string>();
}


inline std::string Settings::getTreeFilename() const
{
    return _parameters.at("treefile").value<std::string>();
}


inline std::string Settings::getSampleProbsFilename() const
{
    return _parameters.at("sampleProbsFilename").value<std::string>();
}


inline double Settings::getGlobalSamplingFraction() const
{
    return _parameters.at("globalSamplingFraction").value<double>();
}


inline double Settings::getUpdateLambdaInitScale() const
{
    return _parameters.at("updateLambdaInitScale").value<double>();
}


inline double Settings::getUpdateMuInitScale() const
{
    return _parameters.at("updateMuInitScale").value<double>();
}


inline double Settings::getUpdateLambdaShiftScale() const
{
    return _parameters.at("updateLambdaShiftScale").value<double>();
}


inline double Settings::getUpdateMuShiftScale() const
{
    return _parameters.at("updateMuShiftScale").value<double>();
}


inline double Settings::getLambdaInit0() const
{
    return _parameters.at("lambdaInit0").value<double>();
}


inline double Settings::getLambdaShift0() const
{
    return _parameters.at("lambdaShift0").value<double>();
}


inline double Settings::getMuInit0() const
{
    return _parameters.at("muInit0").value<double>();
}


inline double Settings::getMuShift0() const
{
    return _parameters.at("muShift0").value<double>();
}


inline double Settings::getUpdateEventRateScale() const
{
    return _parameters.at("updateEventRateScale").value<double>();
}


inline double Settings::getLocalGlobalMoveRatio() const
{
    return _parameters.at("localGlobalMoveRatio").value<double>();
}


inline double Settings::getPoissonRatePrior() const
{
    return _parameters.at("poissonRatePrior").value<double>();
}


inline double Settings::getLambdaInitPrior() const
{
    return _parameters.at("lambdaInitPrior").value<double>();
}


inline double Settings::getLambdaShiftPrior() const
{
    return _parameters.at("lambdaShiftPrior").value<double>();
}


inline double Settings::getMuInitPrior() const
{
    return _parameters.at("muInitPrior").value<double>();
}


inline double Settings::getMuShiftPrior() const
{
    return _parameters.at("muShiftPrior").value<double>();
}


inline double Settings::getLambdaInitRootPrior() const
{
    return _parameters.at("lambdaInitRootPrior").value<double>();
}


inline double Settings::getLambdaShiftRootPrior() const
{
    return _parameters.at("lambdaShiftRootPrior").value<double>();
}


inline double Settings::getMuInitRootPrior() const
{
    return _parameters.at("muInitRootPrior").value<double>();
}


inline double Settings::getMuShiftRootPrior() const
{
    return _parameters.at("muShiftRootPrior").value<double>();
}


inline double Settings::getSegLength() const
{
    return _parameters.at("segLength").value<double>();
}


inline int Settings::getMinCladeSizeForShift() const
{
    return _parameters.at("minCladeSizeForShift").value<int>();
}


inline long int Settings::getSeed() const
{
    return _parameters.at("seed").value<long int>();
}


inline bool Settings::getOverwrite() const
{
    return _parameters.at("overwrite").value<bool>();
}


inline bool Settings::getWriteMeanBranchLengthTrees() const
{
    return _parameters.at("writeMeanBranchLengthTrees").value<bool>();
}


inline std::string Settings::getMCMCoutfile() const
{
    return _parameters.at("mcmcOutfile").value<std::string>();
}


inline std::string Settings::getEventDataOutfile() const
{
    return _parameters.at("eventDataOutfile").value<std::string>();
}


inline std::string Settings::getPriorOutputFileName() const
{
    return _parameters.at("priorOutputFileName").value<std::string>();
}


inline std::string Settings::getLambdaOutfile() const
{
    return _parameters.at("lambdaOutfile").value<std::string>();
}


inline std::string Settings::getMuOutfile() const
{
    return _parameters.at("muOutfile").value<std::string>();
}


inline int Settings::getBranchRatesWriteFreq() const
{
    return _parameters.at("branchRatesWriteFreq").value<int>();
}


inline int Settings::getEventDataWriteFreq() const
{
    return _parameters.at("eventDataWriteFreq").value<int>();
}


inline int Settings::getMCMCwriteFreq() const
{
    return _parameters.at("mcmcWriteFreq").value<int>();
}


inline int Settings::getPrintFreq() const
{
    return _parameters.at("printFreq").value<int>();
}


inline int Settings::getNGENS() const
{
    return _parameters.at("numberGenerations").value<int>();
}


inline double Settings::getUpdateRateEventNumber() const
{
    return _parameters.at("updateRateEventNumber").value<double>();
}


inline double Settings::getUpdateRateEventPosition() const
{
    return _parameters.at("updateRateEventPosition").value<double>();
}


inline double Settings::getUpdateRateEventRate() const
{
    return _parameters.at("updateRateEventRate").value<double>();
}


inline double Settings::getUpdateRateLambda0() const
{
    return _parameters.at("updateRateLambda0").value<double>();
}


inline double Settings::getUpdateRateLambdaShift() const
{
    return _parameters.at("updateRateLambdaShift").value<double>();
}


inline double Settings::getUpdateRateMu0() const
{
    return _parameters.at("updateRateMu0").value<double>();
}


inline double Settings::getUpdateRateMuShift() const
{
    return _parameters.at("updateRateMuShift").value<double>();
}


inline int Settings::getInitialNumberEvents() const
{
    return _parameters.at("initialNumberEvents").value<int>();
}


inline double Settings::getExtinctionProbMax() const
{
    return _parameters.at("ExtinctionProbMax").value<double>();
}


inline std::string Settings::getTraitFile() const
{
    return _parameters.at("traitfile").value<std::string>();
}


inline double Settings::getUpdateBetaInitScale() const
{
    return _parameters.at("updateBetaScale").value<double>();
}


inline double Settings::getUpdateNodeStateScale() const
{
    return _parameters.at("updateNodeStateScale").value<double>();
}


inline double Settings::getUpdateBetaShiftScale() const
{
    return _parameters.at("updateBetaShiftScale").value<double>();
}


inline double Settings::getBetaInit() const
{
    return _parameters.at("betaInit").value<double>();
}


inline double Settings::getBetaShiftInit() const
{
    return _parameters.at("betaShiftInit").value<double>();
}


inline double Settings::getBetaInitPrior() const
{
    return _parameters.at("betaInitPrior").value<double>();
}


inline double Settings::getBetaShiftPrior() const
{
    return _parameters.at("betaShiftPrior").value<double>();
}


inline double Settings::getBetaInitRootPrior() const
{
    return _parameters.at("betaInitRootPrior").value<double>();
}


inline double Settings::getBetaShiftRootPrior() const
{
    return _parameters.at("betaShiftRootPrior").value<double>();
}


inline double Settings::getTraitPriorMin() const
{
    return _parameters.at("traitPriorMin").value<double>();
}


inline double Settings::getTraitPriorMax() const
{
    return _parameters.at("traitPriorMax").value<double>();
}


inline double Settings::getUpdateRateBeta0() const
{
    return _parameters.at("updateRateBeta0").value<double>();
}


inline double Settings::getUpdateRateBetaShift() const
{
    return _parameters.at("updateRateBetaShift").value<double>();
}


inline double Settings::getUpdateRateNodeState() const
{
    return _parameters.at("updateRateNodeState").value<double>();
}


inline bool Settings::getUseObservedMinMaxAsTraitPriors() const
{
    return _parameters.at("useObservedMinMaxAsTraitPriors").value<bool>();
}


inline void Settings::setTraitPriorMin(double x)
{
    _parameters.at("traitPriorMin").setValue<double>(x);
}


inline void Settings::setTraitPriorMax(double x)
{
    _parameters.at("traitPriorMax").setValue<double>(x);
}


inline std::string Settings::getBetaOutfile() const
{
    return _parameters.at("betaOutfile").value<std::string>();
}


inline std::string Settings::getModeltype() const
{
	return _parameters.at("modeltype").value<std::string>();
}


inline void Settings::setUpdateLambdaInitScale(double x)
{
	_parameters.at("updateLambdaInitScale").setValue<double>(x);
}


inline void Settings::setUpdateMuInitScale(double x)
{
	_parameters.at("updateMuInitScale").setValue<double>(x);
}


inline void Settings::setUpdateLambdaShiftScale(double x)
{
	_parameters.at("updateLambdaShiftScale").setValue<double>(x);
}


inline void Settings::setUpdateEventRateScale(double x)
{
	_parameters.at("updateEventRateScale").setValue<double>(x);
}


inline double Settings::getUpdateEventLocationScale() const
{
	return _parameters.at("updateEventLocationScale").value<double>();
}


inline void Settings::setUpdateEventLocationScale(double x)
{
	_parameters.at("updateEventLocationScale").setValue<double>(x);
}


#endif
