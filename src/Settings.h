/*
 *  Settings.h
 *  BAMM
 *
 *  Created by Dan Rabosky on 6/9/12.
 *  Copyright 2012 DLR. All rights reserved.
 *
*/

#ifndef Settings_H
#define Settings_H

#include <string>
#include <stdlib.h>


class Settings
{

private:

    void assertUsingDefault(bool isDefault, const std::string& varName);
    
    bool _allParametersSetToDefaults;

    bool _runTraitModel;
    bool _runSpeciationExtinctionModel;

    bool _sampleFromPriorOnly;
    bool _runMCMC;
    bool _initializeModel;
    bool _loadEventData;
	
	bool _autotune;

    // Files:
    std::string _treefile;
    std::string _eventDataInfile;
	
	std::string _modeltype;

    // Accounting for incomplete sampling:
    bool   _useGlobalSamplingProbability;
    std::string _sampleProbsFilename;
    double _globalSamplingFraction;

    // Parameters relevant to implementation of class Model:
    double _updateLambdaInitScale;
    double _updateMuInitScale;
    double _updateLambdaShiftScale;
    double _updateMuShiftScale;
	double _updateEventLocationScale;
	
    double _lambdaInit0;
    double _lambdaShift0;
    double _muInit0;
    double _muShift0;
	
//   double _MeanSpeciationLengthFraction;
    double _updateEventRateScale;
    double _localGlobalMoveRatio;
    double _poissonRatePrior;

    double _lambdaInitPrior;
    double _lambdaShiftPrior;

    double _muInitPrior;
    double _muShiftPrior;

    double _segLength; // Parm for splitting branches

    int _minCladeSizeForShift;
	
	long int _seed;
	bool _overwrite;
	
    // Class MCMC parameters::General
    std::string _mcmcOutfile;
    std::string _lambdaOutfile;
    std::string _muOutfile;
    std::string _acceptrateOutfile;
    std::string _lambdaNodeOutfile;
    std::string _eventDataOutfile;

    int _treeWriteFreq;
    int _mcmcWriteFreq;
    int _eventDataWriteFreq;
    int _acceptWriteFreq;
    int _printFreq;
    int _NGENS;

    // Class MCMC update weights:
    double _updateRateEventNumber;
    double _updateRateEventPosition;
    double _updateRateEventRate;
    double _updateRateLambda0;
    double _updateRateLambdaShift;
    double _updateRateMu0;
    double _updateRateMuShift;

    // _updateRateNumberTimeVariablePartitions:
    //    if non-zero, does rjMCMC to move between partitions
    //    with fixed and variable rates.
    double _updateRateNumberTimeVariablePartitions;

    // Can fix number of events on tree if _updateRateEventNumber = 0.0
    int _initialNumberEvents;

    // Troubleshooting time-variable model:

    /* Specific parameters for class TraitModel */
    std::string _traitfile;

    double _updateBetaScale;
    double _updateNodeStateScale;
    double _updateBetaShiftScale;

    double _betaInit;
    double _betaShiftInit;

    double _rootPrior;
    double _betaInitPrior;
    double _betaShiftPrior;
    double _traitPriorMin;
    double _traitPriorMax;

    double _updateRateBeta0;
    double _updateRateBetaShift;
    double _updateRateNodeState;

    // Output:
    std::string _betaOutfile;
    std::string _nodeStateOutfile;

    bool _useObservedMinMaxAsTraitPriors;

	// Parameters to create vectors of input variable names 
	//	and associated values: 
	//		these read & populated from controlfile
    std::vector<std::string> _varName;
    std::vector<std::string> _varValue;	
	
    /* ################################################# */
    // Boolean parameters to flag whether default values of a parameter
    // have changed
	
	bool isDefault_modeltype;
	bool isDefault_autotune;
	bool isDefault_treefile;
    bool isDefault_sampleFromPriorOnly;
    bool isDefault_runTraitModel;
    bool isDefault_runSpeciationExtinctionModel;
    bool isDefault_runMCMC;
    bool isDefault_initializeModel;
    bool isDefault_loadEventData;
    bool isDefault_eventDataInfile;
    bool isDefault_sampleProbsFilename;
    bool isDefault_useGlobalSamplingProbability;
    bool isDefault_globalSamplingFraction;
    bool isDefault_updateLambdaInitScale;
    bool isDefault_updateMuInitScale;
    bool isDefault_updateLambdaShiftScale;
    bool isDefault_updateMuShiftScale;
    bool isDefault_lambdaInit0;
    bool isDefault_lambdaShift0;
    bool isDefault_muInit0;
    bool isDefault_muShift0;
    bool isDefault_updateEventRateScale;
    bool isDefault_localGlobalMoveRatio;
    bool isDefault_poissonRatePrior;
    bool isDefault_lambdaInitPrior;
    bool isDefault_lambdaShiftPrior;
    bool isDefault_muInitPrior;
    bool isDefault_muShiftPrior;
    bool isDefault_updateEventLocationScale;
    bool isDefault_segLength;
    bool isDefault_mcmcOutfile;
    bool isDefault_eventDataOutfile;
    bool isDefault_lambdaOutfile;
    bool isDefault_muOutfile;
    bool isDefault_acceptrateOutfile;
    bool isDefault_lambdaNodeOutfile;
    bool isDefault_treeWriteFreq;
    bool isDefault_eventDataWriteFreq;
    bool isDefault_mcmcWriteFreq;
    bool isDefault_acceptWriteFreq;
    bool isDefault_printFreq;
    bool isDefault_NGENS;
    bool isDefault_updateRateEventNumber;
    bool isDefault_updateRateEventPosition;
    bool isDefault_updateRateEventRate;
    bool isDefault_updateRateLambda0;
    bool isDefault_updateRateLambdaShift;
    bool isDefault_updateRateMu0;
    bool isDefault_updateRateMuShift;
    bool isDefault_initialNumberEvents;
    bool isDefault_updateRateNumberTimeVariablePartitions;

    bool isDefault_minCladeSizeForShift;
    
    bool isDefault_seed;
    bool isDefault_overwrite;

    /* specific to trait evolution */
    bool isDefault_traitfile;

    bool isDefault_updateBetaScale;
    bool isDefault_updateNodeStateScale;
    bool isDefault_updateBetaShiftScale;

    bool isDefault_betaInit;
    bool isDefault_betaShift;
    bool isDefault_betaInitPrior;
    bool isDefault_betaShiftPrior;

    bool isDefault_rootPrior;
    bool isDefault_traitPriorMin;
    bool isDefault_traitPriorMax;

    bool isDefault_updateRateBeta0;
    bool isDefault_updateRateNodeState;
    bool isDefault_updateRateBetaShift;

    bool isDefault_useObservedMinMaxAsTraitPriors;
    bool isDefault_betaOutfile;
    bool isDefault_nodeStateOutfile;

public:

    Settings();
    ~Settings();

	void initializeSettingsDevel(std::string controlFilename);
	
	
    void initializeSettingsDefaults_Traits();
    void initializeSettings_Traits();

    void initializeSettings_Diversification();
    void initializeSettingsDefaults_Diversification();

    void checkAreInitialSettingsValid_Diversification();
    void checkAreInitialSettingsValid_Traits();

    void printCurrentSettings_Traits(bool printOnlyChangesToDefaults);
    void printCurrentSettings_Diversification(bool printOnlyChangesToDefaults);

	std::string getModeltype();
	
    bool stringToBool(const char* x);
    void parseCommandLineInput(int argc, std::vector<std::string>& instrings);

    bool areAllParametersSetToDefaults();
    bool getRunTraitModel();
    bool getRunSpeciationExtinctionModel();

    // Functions to access parameters for MCMC/Model/Other
    // from Class Settings object:
    bool getUseGlobalSamplingProbability();
    bool getSampleFromPriorOnly();
    bool getRunMCMC();
    bool getInitializeModel();
	bool getAutotune();

    // Load previous settings?
    bool   getLoadEventData();
    std::string getEventDataInfile();

    // Sampling probabilities:
    std::string getTreeFilename();
    std::string getSampleProbsFilename();
    double getGlobalSamplingFraction();

    // Class Model parameters:
    double getUpdateLambdaInitScale();
    double getUpdateMuInitScale();
    double getUpdateLambdaShiftScale();
    double getUpdateMuShiftScale();
    double getLambdaInit0();
    double getLambdaShift0();
    double getMuInit0();
	double getMuShift0();
	double getUpdateEventLocationScale();

	//    double getMeanSpeciationLengthFraction();
    double getUpdateEventRateScale();
    double getLocalGlobalMoveRatio();
    double getPoissonRatePrior();
    double getLambdaInitPrior();
    double getLambdaShiftPrior();
    double getMuInitPrior();
    double getMuShiftPrior();
    double getSegLength();

    int getMinCladeSizeForShift();
    long int getSeed();
    bool getOverwrite();
	
	// Functions to set MCMC operators for autotuning:
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);
//	void setMeanSpeciationLengthFraction(double x);
	void setUpdateEventRateScale(double x);
	void setUpdateEventLocationScale(double x);
	

    // Class MCMC parameters:
    std::string getMCMCoutfile();
    std::string getEventDataOutfile();
    std::string getLambdaOutfile();
    std::string getMuOutfile();
    std::string getAcceptrateOutfile();
    std::string getLambdaNodeOutfile();

    int getBranchRatesWriteFreq();
    int getEventDataWriteFreq();
    int getMCMCwriteFreq();
    int getAcceptWriteFreq();
    int getPrintFreq();
    int getNGENS();

    // Class MCMC parameter update weights:
    double getUpdateRateEventNumber();
    double getUpdateRateEventPosition();
    double getUpdateRateEventRate();
    double getUpdateRateLambda0();
    double getUpdateRateLambdaShift();
    double getUpdateRateMu0();
    double getUpdateRateMuShift();
    double getUpdateRateNumberTimeVariablePartitions();

    // Other:
    int getInitialNumberEvents();

    /* Parameters specific to trait evolution module */

    std::string getTraitFile();

    double getUpdateBetaScale();
    double getUpdateNodeStateScale();
    double getUpdateBetaShiftScale();
    double getBetaInit();
    double getBetaShiftInit();
    double getRootPrior();
    double getBetaInitPrior();
    double getBetaShiftPrior();
    double getTraitPriorMin();
    double getTraitPriorMax();
    double getUpdateRateBeta0();
    double getUpdateRateBetaShift();
    double getUpdateRateNodeState();
    bool   getUseObservedMinMaxAsTraitPriors();
    void   setTraitPriorMin(double x);
    void   setTraitPriorMax(double x);
    std::string getBetaOutfile();
    std::string getNodeStateOutfile();
};


inline bool Settings::areAllParametersSetToDefaults()
{
    return _allParametersSetToDefaults;
}


inline bool Settings::getRunTraitModel()
{
    return _runTraitModel;
}

inline bool Settings::getAutotune()
{
	return _autotune;
}

inline bool Settings::getRunSpeciationExtinctionModel()
{
    return _runSpeciationExtinctionModel;
}


inline bool Settings::getUseGlobalSamplingProbability()
{
    return _useGlobalSamplingProbability;
}


inline bool Settings::getSampleFromPriorOnly()
{
    return _sampleFromPriorOnly;
}


inline bool Settings::getRunMCMC()
{
    return _runMCMC;
}


inline bool Settings::getInitializeModel()
{
    return _initializeModel;
}


inline bool Settings::getLoadEventData()
{
    return  _loadEventData;
}


inline std::string Settings::getEventDataInfile()
{
    return  _eventDataInfile;
}


inline std::string Settings::getTreeFilename()
{
    return _treefile;
}


inline std::string Settings::getSampleProbsFilename()
{
    return _sampleProbsFilename;
}


inline double Settings::getGlobalSamplingFraction()
{
    return _globalSamplingFraction;
}


inline double Settings::getUpdateLambdaInitScale()
{
    return _updateLambdaInitScale;
}


inline double Settings::getUpdateMuInitScale()
{
    return _updateMuInitScale;
}


inline double Settings::getUpdateLambdaShiftScale()
{
    return _updateLambdaShiftScale;
}


inline double Settings::getUpdateMuShiftScale()
{
    return _updateMuShiftScale;
}


inline double Settings::getLambdaInit0()
{
    return _lambdaInit0;
}


inline double Settings::getLambdaShift0()
{
    return _lambdaShift0;
}


inline double Settings::getMuInit0()
{
    return _muInit0;
}


inline double Settings::getMuShift0()
{
    return _muShift0;
}



inline double Settings::getUpdateEventRateScale()
{
    return _updateEventRateScale;
}


inline double Settings::getLocalGlobalMoveRatio()
{
    return _localGlobalMoveRatio;
}


inline double Settings::getPoissonRatePrior()
{
    return _poissonRatePrior;
}


inline double Settings::getLambdaInitPrior()
{
    return _lambdaInitPrior;
}


inline double Settings::getLambdaShiftPrior()
{
    return _lambdaShiftPrior;
}


inline double Settings::getMuInitPrior()
{
    return _muInitPrior;
}


inline double Settings::getMuShiftPrior()
{
    return _muShiftPrior;
}


inline double Settings::getSegLength()
{
    return _segLength;
}


inline int Settings::getMinCladeSizeForShift()
{
    return _minCladeSizeForShift;
}

inline long int Settings::getSeed()
{
    return _seed;
}

inline bool Settings::getOverwrite()
{
    return _overwrite;
}

inline std::string Settings::getMCMCoutfile()
{
    return _mcmcOutfile;
}


inline std::string Settings::getEventDataOutfile()
{
    return _eventDataOutfile;
}


inline std::string Settings::getLambdaOutfile()
{
    return _lambdaOutfile;
}


inline std::string Settings::getMuOutfile()
{
    return _muOutfile;
}


inline std::string Settings::getAcceptrateOutfile()
{
    return _acceptrateOutfile;
}


inline std::string Settings::getLambdaNodeOutfile()
{
    return _lambdaNodeOutfile;
}


inline int Settings::getBranchRatesWriteFreq()
{
    return _treeWriteFreq;
}


inline int Settings::getEventDataWriteFreq()
{
    return _eventDataWriteFreq;
}


inline int Settings::getMCMCwriteFreq()
{
    return _mcmcWriteFreq;
}


inline int Settings::getAcceptWriteFreq()
{
    return _acceptWriteFreq;
}


inline int Settings::getPrintFreq()
{
    return _printFreq;
}


inline int Settings::getNGENS()
{
    return _NGENS;
}


inline double Settings::getUpdateRateEventNumber()
{
    return _updateRateEventNumber;
}


inline double Settings::getUpdateRateEventPosition()
{
    return _updateRateEventPosition;
}


inline double Settings::getUpdateRateEventRate()
{
    return _updateRateEventRate;
}


inline double Settings::getUpdateRateLambda0()
{
    return _updateRateLambda0;
}


inline double Settings::getUpdateRateLambdaShift()
{
    return _updateRateLambdaShift;
}


inline double Settings::getUpdateRateMu0()
{
    return _updateRateMu0;
}


inline double Settings::getUpdateRateMuShift()
{
    return _updateRateMuShift;
}


inline double Settings::getUpdateRateNumberTimeVariablePartitions()
{
    return _updateRateNumberTimeVariablePartitions;
}


inline int Settings::getInitialNumberEvents()
{
    return _initialNumberEvents;
}


inline std::string Settings::getTraitFile()
{
    return _traitfile;
}


inline double Settings::getUpdateBetaScale()
{
    return _updateBetaScale;
}


inline double Settings::getUpdateNodeStateScale()
{
    return _updateNodeStateScale;
}


inline double Settings::getUpdateBetaShiftScale()
{
    return _updateBetaShiftScale;
}


inline double Settings::getBetaInit()
{
    return _betaInit;
}


inline double Settings::getBetaShiftInit()
{
    return _betaShiftInit;
}


inline double Settings::getRootPrior()
{
    return _rootPrior;
}


inline double Settings::getBetaInitPrior()
{
    return _betaInitPrior;
}


inline double Settings::getBetaShiftPrior()
{
    return _betaShiftPrior;
}


inline double Settings::getTraitPriorMin()
{
    return _traitPriorMin;
}


inline double Settings::getTraitPriorMax()
{
    return _traitPriorMax;
}


inline double Settings::getUpdateRateBeta0()
{
    return _updateRateBeta0;
}


inline double Settings::getUpdateRateBetaShift()
{
    return _updateRateBetaShift;
}


inline double Settings::getUpdateRateNodeState()
{
    return _updateRateNodeState;
}


inline bool Settings::getUseObservedMinMaxAsTraitPriors()
{
    return _useObservedMinMaxAsTraitPriors;
}


inline void Settings::setTraitPriorMin(double x)
{
    _traitPriorMin = x;
}


inline void Settings::setTraitPriorMax(double x)
{
    _traitPriorMax = x;
}


inline std::string Settings::getBetaOutfile()
{
    return _betaOutfile;
}


inline std::string Settings::getNodeStateOutfile()
{
    return _nodeStateOutfile;
}

inline std::string Settings::getModeltype()
{
	return _modeltype;
}

inline void Settings::setUpdateLambdaInitScale(double x)
{
	_updateLambdaInitScale = x;
}

inline void Settings::setUpdateMuInitScale(double x)
{
	_updateMuInitScale = x;
}

inline void Settings::setUpdateLambdaShiftScale(double x)
{
	_updateLambdaShiftScale = x;
}

/* Deprecating
inline void Settings::setMeanSpeciationLengthFraction(double x)
 {
	_MeanSpeciationLengthFraction = x;
}
*/ 
 

inline void Settings::setUpdateEventRateScale(double x){
	_updateEventRateScale = x;
}

inline double Settings::getUpdateEventLocationScale()
{
	return _updateEventLocationScale;
}

inline void Settings::setUpdateEventLocationScale(double x)
{
	_updateEventLocationScale = x;
}

#endif
