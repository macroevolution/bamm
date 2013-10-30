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

using namespace std;



class Settings{

private:
	
	bool			_allParametersSetToDefaults;
	
	bool			_runTraitModel;
	bool			_runSpeciationExtinctionModel;
	
	bool			_sampleFromPriorOnly;
	bool			_runMCMC; 
	bool			_initializeModel;
	bool			_loadEventData;
	
	
	// Files:
	string			_treefile;
	string			_eventDataInfile;
	
	// Accounting for incomplete sampling:
	
	bool			_useGlobalSamplingProbability;	
	string			_sampleProbsFilename;
	double			_globalSamplingFraction;
	
	
	// Parameters relevant to implementation of class Model:
	double			_updateLambdaInitScale;
	double			_updateMuInitScale;
	double			_updateLambdaShiftScale;
	double			_updateMuShiftScale;
	
	double			_lambdaInit0;
	double			_lambdaShift0;
	double			_muInit0;
	double			_muShift0;
	
	double			_MeanSpeciationLengthFraction;
	double			_updateEventRateScale;
	double			_localGlobalMoveRatio;
	double			_targetNumber;
	
	double			_lambdaInitPrior;
	double			_lambdaShiftPrior;
	
	double			_muInitPrior;
	double			_muShiftPrior;
	
	double			_segLength; // Parm for splitting branches...
	
	int				_minCladeSizeForShift;
	
	// Class MCMC parameters::General
	string			_mcmcOutfile;
	string			_lambdaOutfile;
	string			_muOutfile;
	string			_acceptrateOutfile;
	string			_lambdaNodeOutfile;
	string			_eventDataOutfile;
	
	
	
	int				_treeWriteFreq;
	int				_mcmcWriteFreq;
	int				_eventDataWriteFreq;
	int				_acceptWriteFreq;
	int				_printFreq;
	int				_NGENS;
	
	// Class MCMC update weights:
	double			_updateRateEventNumber;
	double			_updateRateEventPosition;
	double			_updateRateEventRate;
	double			_updateRateLambda0;
	double			_updateRateLambdaShift;
	double			_updateRateMu0;
	double			_updateRateMuShift;
	
	// _updateRateNumberTimeVariablePartitions:
	//		if non-zero, does rjMCMC to move between partitions 
	//		with fixed and variable rates.
	double			_updateRateNumberTimeVariablePartitions;

	// Other:
	int				_initialNumberEvents; // can fix number of events on tree if _updateRateEventNumber = 0.0

	// Troubleshooting time-variable model:

	/* Specific parameters for class TraitModel */
	string			_traitfile;
	
	double			_updateBetaScale;
	double			_updateNodeStateScale;
	double			_updateBetaShiftScale;
	
	double			_betaInit;
	double			_betaShiftInit;
	
	double			_rootPrior;
	double			_betaInitPrior;
	double			_betaShiftPrior;
	double			_traitPriorMin;
	double			_traitPriorMax;
	
	double			_updateRateBeta0;
	double			_updateRateBetaShift;
	double			_updateRateNodeState;
	
	// Output:
	string			_betaOutfile;
	string			_nodeStateOutfile;
	
	bool			_useObservedMinMaxAsTraitPriors;
	
	
	/* ################################################# */ 
	// Boolean parameters to flag whether default values of a parameter have changed
	bool	isDefault_treefile;
	bool	isDefault_sampleFromPriorOnly;
	bool	isDefault_runTraitModel;
	bool	isDefault_runSpeciationExtinctionModel;
	bool	isDefault_runMCMC;
	bool	isDefault_initializeModel;
	bool	isDefault_loadEventData;
	bool	isDefault_eventDataInfile;
	bool	isDefault_sampleProbsFilename;
	bool	isDefault_useGlobalSamplingProbability;
	bool	isDefault_globalSamplingFraction;
	bool	isDefault_updateLambdaInitScale;
	bool	isDefault_updateMuInitScale;
	bool	isDefault_updateLambdaShiftScale;
	bool	isDefault_updateMuShiftScale;
	bool	isDefault_lambdaInit0;
	bool	isDefault_lambdaShift0;
	bool	isDefault_muInit0;
	bool	isDefault_muShift0;
	bool	isDefault_updateEventRateScale;
	bool	isDefault_localGlobalMoveRatio;
	bool	isDefault_targetNumber;
	bool	isDefault_lambdaInitPrior;
	bool	isDefault_lambdaShiftPrior;
	bool	isDefault_muInitPrior;
	bool	isDefault_muShiftPrior;
	bool	isDefault_MeanSpeciationLengthFraction;
	bool	isDefault_segLength;
	bool	isDefault_mcmcOutfile;
	bool	isDefault_eventDataOutfile;
	bool	isDefault_lambdaOutfile;
	bool	isDefault_muOutfile;
	bool	isDefault_acceptrateOutfile;
	bool	isDefault_lambdaNodeOutfile;
	bool	isDefault_treeWriteFreq;
	bool	isDefault_eventDataWriteFreq;
	bool	isDefault_mcmcWriteFreq;
	bool	isDefault_acceptWriteFreq;
	bool	isDefault_printFreq;
	bool	isDefault_NGENS;
	bool	isDefault_updateRateEventNumber;
	bool	isDefault_updateRateEventPosition;
	bool	isDefault_updateRateEventRate;
	bool	isDefault_updateRateLambda0;
	bool	isDefault_updateRateLambdaShift;
	bool	isDefault_updateRateMu0;
	bool	isDefault_updateRateMuShift;
	bool	isDefault_initialNumberEvents;	
	bool	isDefault_updateRateNumberTimeVariablePartitions;
	
	bool	isDefault_minCladeSizeForShift;
	
	/* specific to trait evolution */
	bool	isDefault_traitfile;
	
	bool	isDefault_updateBetaScale;
	bool	isDefault_updateNodeStateScale;
	bool	isDefault_updateBetaShiftScale;
	
	bool	isDefault_betaInit;
	bool	isDefault_betaShift;
	bool	isDefault_betaInitPrior;
	bool	isDefault_betaShiftPrior;
	
	bool	isDefault_rootPrior;
	bool	isDefault_traitPriorMin;
	bool	isDefault_traitPriorMax;
	
	bool	isDefault_updateRateBeta0;
	bool	isDefault_updateRateNodeState;
	bool	isDefault_updateRateBetaShift;
	
	bool	isDefault_useObservedMinMaxAsTraitPriors;
	bool	isDefault_betaOutfile;
	bool	isDefault_nodeStateOutfile;
	
	
public:
	Settings(void);
	~Settings(void);

	void	trait_initializeSettings(void);
	void	trait_initializeSettings(string controlFilename);
	
	void	initializeSettings(string controlFilename);
	
	void	initializeSettings(void);
	
	void	trait_printCurrentSettings(bool printOnlyChangesToDefaults);
	void	printCurrentSettings(bool printOnlyChangesToDefaults);

	

	bool	stringToBool(const char * x);
	void	parseCommandLineInput(int argc, vector<string> &instrings, string modeltype);
	
	bool	areAllParametersSetToDefaults(void)		{ return _allParametersSetToDefaults;		}
	bool	getRunTraitModel(void)					{ return _runTraitModel;					}
	bool	getRunSpeciationExtinctionModel(void)	{ return _runSpeciationExtinctionModel;		}
		
	
// Functions to access parameters for MCMC/Model/Other from Class Settings object:	
	bool	getUseGlobalSamplingProbability(void)	{  return _useGlobalSamplingProbability;	}
	bool	getSampleFromPriorOnly(void)			{  return _sampleFromPriorOnly;				}
	bool	getRunMCMC(void)						{  return _runMCMC;							}
	bool	getInitializeModel(void)				{  return _initializeModel;					}
	
// Load previous settings?
	bool	getLoadEventData(void)					{ return  _loadEventData;					}
	string	getEventDataInfile(void)				{ return  _eventDataInfile;					}
	
// Sampling probabilities:

	string	getTreeFilename(void)					{  return _treefile;						}
	string	getSampleProbsFilename(void)			{  return _sampleProbsFilename;				}
	double	getGlobalSamplingFraction(void)			{  return _globalSamplingFraction;			}
	
	// Class Model parameters:
	double	getUpdateLambdaInitScale(void)			{  return _updateLambdaInitScale;			}
	double	getUpdateMuInitScale(void)				{  return _updateMuInitScale;				}
	double	getUpdateLambdaShiftScale(void)			{  return _updateLambdaShiftScale;			}
	double	getUpdateMuShiftScale(void)				{  return _updateMuShiftScale;				}
	double	getLambdaInit0(void)					{  return _lambdaInit0;						}
	double	getLambdaShift0(void)					{  return _lambdaShift0;					}
	double	getMuInit0(void)						{  return _muInit0;							}
	double	getMuShift0(void)						{  return _muShift0;						}
	double	getMeanSpeciationLengthFraction(void)	{  return _MeanSpeciationLengthFraction;	}
	double	getUpdateEventRateScale(void)			{  return _updateEventRateScale;			}
	double	getLocalGlobalMoveRatio(void)			{  return _localGlobalMoveRatio;			}
	double	getTargetNumberOfEvents(void)			{  return _targetNumber;					}
	double	getLambdaInitPrior(void)				{  return _lambdaInitPrior;					}
	double	getLambdaShiftPrior(void)				{  return _lambdaShiftPrior;				}
	double	getMuInitPrior(void)					{  return _muInitPrior;						}
	double	getMuShiftPrior(void)					{  return _muShiftPrior;					}
	double	getSegLength(void)						{  return _segLength;						}
	
	int		getMinCladeSizeForShift(void)			{  return _minCladeSizeForShift;			}
	
	
	// Class MCMC parameters:
	string	getMCMCoutfile(void)					{  return _mcmcOutfile;						}
	string	getEventDataOutfile(void)				{  return _eventDataOutfile;				}
	string	getLambdaOutfile(void)					{  return _lambdaOutfile;					}
	string	getMuOutfile(void)						{  return _muOutfile;						}
	string	getAcceptrateOutfile(void)				{  return _acceptrateOutfile;				}
	string  getLambdaNodeOutfile(void)				{  return _lambdaNodeOutfile;				}
	
	int		getTreeWriteFreq(void)					{  return _treeWriteFreq;					}
	int		getEventDataWriteFreq(void)				{  return _eventDataWriteFreq;				}
	int		getMCMCwriteFreq(void)					{  return _mcmcWriteFreq;					}
	int		getAcceptWriteFreq(void)				{  return _acceptWriteFreq;					}
	int		getPrintFreq(void)						{  return _printFreq;						}
	int		getNGENS(void)							{  return _NGENS;							}
	
	
	// Class MCMC parameter update weights:
	double	getUpdateRateEventNumber(void)			{  return _updateRateEventNumber;			}
	double	getUpdateRateEventPosition(void)		{  return _updateRateEventPosition;			}
	double	getUpdateRateEventRate(void)			{  return _updateRateEventRate;				}
	double	getUpdateRateLambda0(void)				{  return _updateRateLambda0;				}
	double	getUpdateRateLambdaShift(void)			{  return _updateRateLambdaShift;			}
	double	getUpdateRateMu0(void)					{  return _updateRateMu0;					}
	double	getUpdateRateMuShift(void)				{  return _updateRateMuShift;				}
	double	getUpdateRateNumberTimeVariablePartitions(void)	{ return _updateRateNumberTimeVariablePartitions;	}
	
	
	// Other:
	int		getInitialNumberEvents(void)			{  return _initialNumberEvents;				}

	/* Parameters specific to trait evolution module */ 
	
	string			getTraitFile(void)				{ return _traitfile;						}
	
	double			getUpdateBetaScale(void)		{ return _updateBetaScale;					}
	double			getUpdateNodeStateScale(void)	{ return _updateNodeStateScale;				}
	double			getUpdateBetaShiftScale(void)	{ return _updateBetaShiftScale;				}
	
	double			getBetaInit(void)				{ return _betaInit;							}
	double			getBetaShiftInit(void)			{ return _betaShiftInit;					}
	
	double			getRootPrior(void)				{ return _rootPrior;						}
	double			getBetaInitPrior(void)			{ return _betaInitPrior;					} 
	double			getBetaShiftPrior(void)			{ return _betaShiftPrior;					}
	double			getTraitPriorMin(void)			{ return _traitPriorMin;					}
	double			getTraitPriorMax(void)			{ return _traitPriorMax;					}
	
	double			getUpdateRateBeta0(void)		{ return _updateRateBeta0;					}
	double			getUpdateRateBetaShift(void)	{ return _updateRateBetaShift;				}
	double			getUpdateRateNodeState(void)	{ return _updateRateNodeState;				}
	
	bool			getUseObservedMinMaxAsTraitPriors(void)	{ return _useObservedMinMaxAsTraitPriors;	}
	
	void			setTraitPriorMin(double x)		{ _traitPriorMin = x;			}
	void			setTraitPriorMax(double x)		{ _traitPriorMax = x;			}
	
	string			getBetaOutfile(void)			{ return _betaOutfile;			}
	string			getNodeStateOutfile(void)		{ return _nodeStateOutfile;		}
	
};





#endif
