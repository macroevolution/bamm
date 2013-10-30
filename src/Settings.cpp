/*
 *  settings.cpp
 *  rateshiftcombined
 *
 *  Created by Dan Rabosky on 2/6/12.
  *
 */

 
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cctype>

#include "Settings.h"

using namespace std;

/*
	Should rewrite this to NOT initialize any variables - just set to 0.
 
 */

Settings::Settings(void){
	
	_allParametersSetToDefaults = false;
	// Parameters used in main()
	_treefile = "EMPTY_STRING";
	_sampleProbsFilename = "EMPTY_STRING";
	_eventDataInfile = "EMPTY_STRING";
	
	_runTraitModel = false;
	_runSpeciationExtinctionModel = false;
	_sampleFromPriorOnly = false;
	_runMCMC = false;
	_initializeModel = false;
	_loadEventData = false;
	
	_useGlobalSamplingProbability = true;
	_globalSamplingFraction = NULL;
	
	// Class Model parameters:
	_updateLambdaInitScale = NULL;
	_updateMuInitScale = NULL;
	_updateLambdaShiftScale = NULL;
	_updateMuShiftScale = NULL;
	_lambdaInit0 = NULL;
	_lambdaShift0 = NULL;
	_muInit0 = NULL;
	_muShift0 = NULL;
	_updateEventRateScale = NULL;
	_localGlobalMoveRatio = NULL;
	_targetNumber = NULL;
	_lambdaInitPrior = NULL;
	_lambdaShiftPrior = NULL;
	_muInitPrior = NULL;
	_muShiftPrior = 1.0;  // This is only set for convenience for now:
	_MeanSpeciationLengthFraction = NULL; 
	_segLength = NULL;
	
	_minCladeSizeForShift = 1;
	
	// Parameters for implementation of class MCMC:
	_mcmcOutfile			=		"BAMM_mcmc_out.txt";
	_eventDataOutfile		=		"BAMM_eventdata.txt";
	_lambdaOutfile			=		"BAMM_lambda_rates.txt";
	_muOutfile				=		"BAMM_mu_rates.txt";
	_acceptrateOutfile		=		"BAMM_mcmc_accept.txt";
	_lambdaNodeOutfile		=		"BAMM_nodeLambda.txt";
	
	_treeWriteFreq			=		NULL;
	_eventDataWriteFreq		=		NULL;
	_mcmcWriteFreq			=		NULL;
	_acceptWriteFreq		=		NULL;
	_printFreq				=		NULL;
	_NGENS					=		NULL;
	
	_updateRateEventNumber					=	NULL;
	_updateRateEventPosition				=	NULL;
	_updateRateEventRate					=	NULL;
	_updateRateLambda0						=	NULL;
	_updateRateLambdaShift					=	NULL;
	_updateRateMu0							=	NULL;
	_updateRateMuShift						=	NULL;
	_updateRateNumberTimeVariablePartitions	=	NULL;
	
	// Other:
	_initialNumberEvents					=	0;	

	/*******************************************/
	// Set Boolean indicator variables to track
	//	 changes to default parameters:
	
	isDefault_treefile									= true;
	isDefault_runTraitModel								= true;
	isDefault_runSpeciationExtinctionModel				= true;
	isDefault_sampleFromPriorOnly						= true;
	isDefault_runMCMC									= true;
	isDefault_initializeModel							= true;
	isDefault_loadEventData								= true;
	isDefault_useGlobalSamplingProbability				= true;
	isDefault_sampleProbsFilename						= true;
	isDefault_globalSamplingFraction					= true;
	isDefault_updateLambdaInitScale						= true;
	isDefault_updateMuInitScale							= true;
	isDefault_updateLambdaShiftScale					= true;
	isDefault_updateMuShiftScale						= true;
	isDefault_lambdaInit0								= true;
	isDefault_lambdaShift0								= true;
	isDefault_muInit0									= true;
	isDefault_muShift0									= true;
	isDefault_updateEventRateScale						= true;
	isDefault_localGlobalMoveRatio						= true;
	isDefault_targetNumber								= true;
	isDefault_lambdaInitPrior							= true;
	isDefault_lambdaShiftPrior							= true;
	isDefault_muInitPrior								= true;
	isDefault_muShiftPrior								= true;
	isDefault_MeanSpeciationLengthFraction				= true;
	isDefault_segLength									= true;
	isDefault_mcmcOutfile								= true;
	isDefault_eventDataOutfile							= true;
	isDefault_lambdaOutfile								= true;
	isDefault_muOutfile									= true;
	isDefault_acceptrateOutfile							= true;
	isDefault_lambdaNodeOutfile							= true;
	isDefault_treeWriteFreq								= true;
	isDefault_eventDataWriteFreq						= true;
	isDefault_mcmcWriteFreq								= true;
	isDefault_acceptWriteFreq							= true;
	isDefault_printFreq									= true;
	isDefault_NGENS										= true;
	isDefault_updateRateEventNumber						= true;
	isDefault_updateRateEventPosition					= true;
	isDefault_updateRateEventRate						= true;
	isDefault_updateRateLambda0							= true;
	isDefault_updateRateLambdaShift						= true;
	isDefault_updateRateMu0								= true;
	isDefault_updateRateMuShift							= true;
	isDefault_updateRateNumberTimeVariablePartitions	= true;
	isDefault_initialNumberEvents						= true;
	
	isDefault_eventDataInfile							= true;
	isDefault_minCladeSizeForShift						= true;
	
	// End block of Booleans
	/*******************************************/
	
	/* Trait evolution module parameters  */
	_traitfile					=	"EMPTY";
	_updateBetaScale			=	NULL;
	_updateNodeStateScale		=	NULL;
	_updateBetaShiftScale		=	NULL;
	
	_betaInit					=	NULL;
	_betaShiftInit				=	NULL;
	
	_rootPrior					=	NULL;
	_betaInitPrior				=	NULL;
	_betaShiftPrior				=	NULL;
	
	_traitPriorMin				=	NULL;
	_traitPriorMax				=	NULL;
	
	isDefault_traitfile			= true;
	isDefault_updateBetaScale	= true;
	isDefault_updateNodeStateScale = true;
	isDefault_betaInit			= true;
	isDefault_rootPrior			= true;
	isDefault_traitPriorMin		= true;
	isDefault_traitPriorMax		= true;	
	isDefault_useObservedMinMaxAsTraitPriors = true;
	isDefault_betaOutfile		= true;
	isDefault_nodeStateOutfile	= true;

	isDefault_updateRateBetaShift	= true;
	isDefault_updateRateBeta0		= true;
	isDefault_updateRateNodeState	= true;
	
	_useObservedMinMaxAsTraitPriors	= true;
	
	_betaOutfile		=	"BAMM_beta_rates.txt";
	_nodeStateOutfile	=	"BAMM_nodestates.txt";
	
}


Settings::~Settings(void){
	
 


}




/*
	Initializes all settings directly in the compiled code.

 
 */

void Settings::initializeSettings(void){

	_allParametersSetToDefaults = true;
	
	_runSpeciationExtinctionModel	= true;
	_runTraitModel					= false;
	
	// Parameters used in main()
	_treefile = "test_tree.tre";

	_sampleFromPriorOnly = false;
	_runMCMC = true;
	_initializeModel = true;
	
	_useGlobalSamplingProbability	= true;
	_sampleProbsFilename = "skinksprobs.txt";
	
	_globalSamplingFraction			= 0.95;
	
	// Class Model parameters:
	_updateLambdaInitScale			= 2.0;
	_updateMuInitScale				= 2.0;
	_updateLambdaShiftScale			= 2.0;
	_updateMuShiftScale				= 0.0; // 0
	_lambdaInit0					= 1.0;
	_lambdaShift0					= 0.0;
	_muInit0						= 0.1;
	_muShift0						= 0.0;
	_updateEventRateScale			= 4.0;
	_localGlobalMoveRatio			= 10.0;
	_targetNumber					= 1.0;
	_lambdaInitPrior				= 1.0;
	_lambdaShiftPrior				= 0.5;
	_muInitPrior					= 1.0;
	_muShiftPrior					= 0.5; // 0
	_MeanSpeciationLengthFraction	= 0.2; 
	_segLength						= 1.0;
	
	_minCladeSizeForShift			= 1;
	
	// Parameters for implementation of class MCMC:
	_mcmcOutfile			=		"mcmc_out.txt";
	_eventDataOutfile		=		"eventdata.txt";
	_eventDataInfile		=		"EMPTY_STRING";
	_lambdaOutfile			=		"lambda_rates.txt";
	_muOutfile				=		"mu_rates.txt";
	_acceptrateOutfile		=		"mcmc_accept.txt";
	_lambdaNodeOutfile		=		"nodeLambda.txt";
	
	
	_treeWriteFreq			=		5000;
	_eventDataWriteFreq		=		5000;
	_mcmcWriteFreq			=		1000;
	_acceptWriteFreq		=		1000;
	_printFreq				=		1000;
	_NGENS					=		2000000;
	
	_updateRateEventNumber		=	1.0;
	_updateRateEventPosition	=	1.0;
	_updateRateEventRate		=	1.0;
	_updateRateLambda0			=	10.0;
	_updateRateLambdaShift		=	10.0;
	_updateRateMu0				=	10.0;
	_updateRateMuShift			=	0.0; // 0.0
	
	_updateRateNumberTimeVariablePartitions	= 0.0;
	
	// Other:
	_initialNumberEvents = 0;
	

}

void Settings::trait_initializeSettings(void){
	_allParametersSetToDefaults = true;
	
	_runSpeciationExtinctionModel	= false;
	_runTraitModel					= true;
	
	// Parameters used in main()
	_treefile = "test_tree.txt";
	_traitfile	= "morph.txt";
 	
	_sampleFromPriorOnly = false;
	_runMCMC = true;
	_initializeModel = true;
 
	// Class Model parameters:
	_updateBetaScale				= 0.25;
	_updateNodeStateScale			= 0.25;
	_updateBetaShiftScale			= 0.25;
	
	_updateEventRateScale			= 2.0;
	_localGlobalMoveRatio			= 10.0;
	_targetNumber					= 0.5;
	_MeanSpeciationLengthFraction	= 0.2; 
	
	_betaInit			=	0.1;
	_betaShiftInit		=	0.0;
	
	_betaInitPrior		=	1.0;
	_betaShiftPrior		=	0.1;
	
	// These still have default values: 
	_useObservedMinMaxAsTraitPriors	= true;
	_traitPriorMin					=	NULL;
	_traitPriorMax					=	NULL;	
	
	// Parameters for implementation of class MCMC:
	_mcmcOutfile			=		"BAMMt_mcmc_out.txt";
	_eventDataOutfile		=		"BAMMt_eventdata.txt";
	_betaOutfile			=		"BAMMt_beta_rates.txt";
	_nodeStateOutfile		=		"BAMMt_nodestates.txt";
	_acceptrateOutfile		=		"BAMMt_mcmc_accept.txt";
	
	_treeWriteFreq			=		50000;
	_eventDataWriteFreq		=		50000;
	_mcmcWriteFreq			=		1000;
	_acceptWriteFreq		=		1000;
	_printFreq				=		10000;
	_NGENS					=		10000000;
	
	_updateRateEventNumber		=	1.0;
	_updateRateEventPosition	=	1.0;
	_updateRateEventRate		=	1.0;
	
	_updateRateBeta0			=	1.0;
	_updateRateNodeState		=	25.0;
	_updateRateBetaShift		=	1.0;
	
	_updateRateNumberTimeVariablePartitions	= 0.0;
	
	// Other:
	_initialNumberEvents = 0;
	
	
}


// Fixing this : March 4 2013
void Settings::trait_initializeSettings(string controlFilename){
	
	
	
	_allParametersSetToDefaults = false;
	
	int parmCount = 0;
	int ppw = 35;
	
	ifstream infile(controlFilename.c_str());
	if (!infile.good()){
		cout << "Control filename invalid" << endl;
		cout << "Exiting." << endl;
		throw;
	}
	
	cout << "Reading control file <<" << controlFilename.c_str() << ">>" << endl;
	
	vector<string> varName;
	vector<string> varValue;
	vector<string> paramsNotFound;
	string s1, s2;
	
	
	while (infile){
		getline(infile, s1, '\n');
		
		// strip whitespace out of tempstring:
		//		both spaces and tabs:
		
		// What is the int(*)(int) doing????
		s1.erase(remove_if(s1.begin(), s1.end(), (int(*)(int))isspace), s1.end());
		
		//tempstring.erase(remove_if(tempstring.begin(), tempstring.end(), isspace), tempstring.end());
		
		
		// Only add if has size > 0 (gets rid of empty lines)
		if (s1.size() > 0){
			vector<string> tmpstr;
			
			// NOw use second getline to split by '=' characters:
			
			istringstream stemp(s1);
			while (getline(stemp, s2, '=')){
				tmpstr.push_back(s2);
			}
			if (tmpstr.size() == 2){
				varName.push_back(tmpstr[0]);
				varValue.push_back(tmpstr[1]);
			}else{
				cout << "Invalid size of input line in control file" << endl;
			}
			
			//getline(stemp, s2, '=');
			//cout << s2 << endl;
			//stringvec.push_back(s1);
			
		}
		
		if (infile.peek() == EOF)
			break;
	}
	
	
	for (int i = 0; i < varName.size(); i++){
		//cout << setw(30) << varName[i] << setw(20) << varValue[i] << endl;
		
		if (varName[i] == "treefile"){
			_treefile = varValue[i];
			isDefault_treefile = false;
 		}else if (varName[i] == "traitfile"){
			_traitfile = varValue[i];
			isDefault_traitfile = false;
		}else if(varName[i] == "runSpeciationExtinctionModel"){
			
		}else if(varName[i] == "runTraitModel"){
			
		}else if (varName[i] == "sampleFromPriorOnly"){
			_sampleFromPriorOnly = stringToBool(varValue[i].c_str());
			
			isDefault_sampleFromPriorOnly = false;
		}else if (varName[i] == "initializeModel"){
			_initializeModel = stringToBool(varValue[i].c_str());
			isDefault_initializeModel = false;
		}else if (varName[i] == "runMCMC"){
			_runMCMC = stringToBool(varValue[i].c_str());
			isDefault_runMCMC = false;
		}else if (varName[i] == "updateBetaScale"){
			_updateBetaScale = atof(varValue[i].c_str());
			isDefault_updateBetaScale = false;
			//cout << left << setw(ppw) << "updateLambdaInitScale" << "\t" << _updateLambdaInitScale << endl;
			
		}else if (varName[i] == "updateNodeStateScale"){
			_updateNodeStateScale = atof(varValue[i].c_str());
			isDefault_updateNodeStateScale = false;
		}else if (varName[i] == "updateBetaShiftScale"){
			_updateBetaShiftScale = atof(varValue[i].c_str());
			isDefault_updateBetaShiftScale = false;
		}else if (varName[i] == "betaInit"){
			_betaInit = atof(varValue[i].c_str());
			isDefault_betaInit = false;
		}else if (varName[i] == "betaShiftInit"){
			_betaShiftInit = atof(varValue[i].c_str());
			isDefault_betaShift = false;			 
		}else if (varName[i] == "updateEventRateScale"){
			_updateEventRateScale = atof(varValue[i].c_str());
			isDefault_updateEventRateScale = false;
		}else if (varName[i] == "localGlobalMoveRatio"){
			_localGlobalMoveRatio = atof(varValue[i].c_str());	
			isDefault_localGlobalMoveRatio = false;
		}else if (varName[i] == "targetNumber"){
			_targetNumber = atof(varValue[i].c_str());
			isDefault_targetNumber = false;
		}else if (varName[i] == "betaInitPrior"){
			_betaInitPrior = atof(varValue[i].c_str());
			isDefault_betaInitPrior = false;
		}else if (varName[i] == "betaShiftPrior"){
			_betaShiftPrior = atof(varValue[i].c_str());
			isDefault_betaShiftPrior = false;
		}else if (varName[i] == "useObservedMinMaxAsTraitPriors"){
			_useObservedMinMaxAsTraitPriors = stringToBool(varValue[i].c_str());
			isDefault_useObservedMinMaxAsTraitPriors = false;    
		}else if (varName[i] == "traitPriorMin"){
			_traitPriorMin = atof(varValue[i].c_str());
			isDefault_traitPriorMin = false;
		}else if (varName[i] == "traitPriorMax"){
			_traitPriorMax = atof(varValue[i].c_str());
			isDefault_traitPriorMax = false;
		}else if (varName[i] == "mcmcOutfile"){
			_mcmcOutfile = varValue[i];
			isDefault_mcmcOutfile = false;
		}else if (varName[i] == "eventDataOutfile"){
			_eventDataOutfile = varValue[i];
			isDefault_eventDataOutfile = false;
		}else if (varName[i] == "betaOutfile"){
			_betaOutfile = varValue[i];
			isDefault_betaOutfile = false;
		}else if (varName[i] == "nodeStateOutfile"){
			_nodeStateOutfile = varValue[i];
			isDefault_nodeStateOutfile = false;
		}else if (varName[i] == "acceptrateOutfile"){
			_acceptrateOutfile = varValue[i];
			isDefault_acceptrateOutfile = false;
		}else if (varName[i] == "treeWriteFreq"){
			_treeWriteFreq = atoi(varValue[i].c_str());
			isDefault_treeWriteFreq = false;
		}else if (varName[i] == "eventDataWriteFreq"){
			_eventDataWriteFreq = atoi(varValue[i].c_str());
			isDefault_eventDataWriteFreq = false;
		}else if (varName[i] == "mcmcWriteFreq"){
			_mcmcWriteFreq = atoi(varValue[i].c_str());
			isDefault_mcmcWriteFreq = false;
		}else if (varName[i] == "acceptWriteFreq"){
			_acceptWriteFreq = atoi(varValue[i].c_str());
			isDefault_acceptWriteFreq = false;
		}else if (varName[i] == "printFreq"){
			_printFreq = atoi(varValue[i].c_str());
			isDefault_printFreq = false;
		}else if (varName[i] == "NumberGenerations"){
			_NGENS = atoi(varValue[i].c_str());
			isDefault_NGENS = false;
		}else if (varName[i] == "updateRateEventNumber"){
			_updateRateEventNumber  = atof(varValue[i].c_str());
			isDefault_updateRateEventNumber = false;
		}else if (varName[i] == "updateRateEventPosition"){
			_updateRateEventPosition  = atof(varValue[i].c_str());
			isDefault_updateRateEventPosition = false;
		}else if (varName[i] == "updateRateEventRate"){
			_updateRateEventRate = atof(varValue[i].c_str());
			isDefault_updateRateEventRate = false;
		}else if (varName[i] == "updateRateBeta0"){
			_updateRateBeta0 = atof(varValue[i].c_str());
			isDefault_updateRateBeta0 = false;
		}else if (varName[i] == "updateRateBetaShift"){
			_updateRateBetaShift  = atof(varValue[i].c_str());
			isDefault_updateRateBetaShift = false;
		}else if (varName[i] == "updateRateNodeState"){
			_updateRateNodeState = atoi(varValue[i].c_str());
			isDefault_updateRateNodeState = false;
		}else if (varName[i] == "initialNumberEvents"){
			_initialNumberEvents  = atoi(varValue[i].c_str());
			isDefault_initialNumberEvents = false;
		}else if(varName[i] == "loadEventData" ){
			_loadEventData = stringToBool(varValue[i].c_str());	
			isDefault_loadEventData = false;
		}else if (varName[i] == "eventDataInfile"){
			_eventDataInfile = varValue[i].c_str();
			isDefault_eventDataInfile = false;
		}else{
			// Parameter not found:
			//		add to list of potentially bad/misspelled params 
			//		and print for user.
			paramsNotFound.push_back(varName[i]);
		}	
		
	}
	
	cout << "Read a total of <<" << varName.size() << ">> parameter settings from control file" << endl;
	if (paramsNotFound.size() > 0){
		cout << endl << "********************************" << endl;
		cout << "BAMM error: one or more parameters from control file do not correspond" << endl;
		cout << "\tto valid model parameters.";
		cout << "Check the following to see if they are ";
		cout << "\tspecified (or spelled) correctly:" << endl << endl;
		for (int i = 0; i < paramsNotFound.size(); i++){
			cout << setw(30) << paramsNotFound[i] << endl;		
		}
		cout << endl << "********************************" << endl << endl;
		cout << "Execution of BAMM terminated..." << endl;
		exit(1);
		
	}
	
	//	Here we have a print block to output Settings:
	//	Any parameters NOT set will have the defaults.
	//	Thus user can specify a control file with ONLY 
	//		those parameters that they wish to change 
	//		from the defaults, eg inputfilename etc.
	
	//	Output list of default parameters.
	

}

void Settings::trait_printCurrentSettings(bool printOnlyChangesToDefaults){
	cout << "print settings for trait module not yet supported" << endl;
	exit(1);
}

void Settings::initializeSettings(string controlFilename){

	_allParametersSetToDefaults = false;
	
	int parmCount = 0;
	int ppw = 35;
	
	ifstream infile(controlFilename.c_str());
	if (!infile.good()){
		cout << "Control filename invalid" << endl;
		cout << "Exiting." << endl;
		throw;
	}
	
	cout << "Reading control file <<" << controlFilename.c_str() << ">>" << endl;
	
	vector<string> varName;
	vector<string> varValue;
	vector<string> paramsNotFound;
	string s1, s2;
		
	
	while (infile){
		getline(infile, s1, '\n');
		
		// strip whitespace out of tempstring:
		//		both spaces and tabs:
 
		// What is the int(*)(int) doing????
		s1.erase(remove_if(s1.begin(), s1.end(), (int(*)(int))isspace), s1.end());
 
		//tempstring.erase(remove_if(tempstring.begin(), tempstring.end(), isspace), tempstring.end());
		
		
		// Only add if has size > 0 (gets rid of empty lines)
		if (s1.size() > 0){
			vector<string> tmpstr;
			
			// NOw use second getline to split by '=' characters:
			
			istringstream stemp(s1);
			while (getline(stemp, s2, '=')){
				tmpstr.push_back(s2);
			}
			if (tmpstr.size() == 2){
				varName.push_back(tmpstr[0]);
				varValue.push_back(tmpstr[1]);
			}else{
				cout << "Invalid size of input line in control file" << endl;
			}
			
			//getline(stemp, s2, '=');
			//cout << s2 << endl;
			//stringvec.push_back(s1);
		
		}
	
		if (infile.peek() == EOF)
			break;
	}
 

	for (int i = 0; i < varName.size(); i++){
		//cout << setw(30) << varName[i] << setw(20) << varValue[i] << endl;
 
		if (varName[i] == "treefile"){
			_treefile = varValue[i];
			isDefault_treefile = false;
 		}else if(varName[i] == "runSpeciationExtinctionModel"){
		
		}else if(varName[i] == "runTraitModel"){
		
		}else if (varName[i] == "sampleFromPriorOnly"){
			_sampleFromPriorOnly = stringToBool(varValue[i].c_str());
			
			isDefault_sampleFromPriorOnly = false;
		}else if (varName[i] == "initializeModel"){
			_initializeModel = stringToBool(varValue[i].c_str());
			isDefault_initializeModel = false;
		}else if (varName[i] == "runMCMC"){
			_runMCMC = stringToBool(varValue[i].c_str());
			isDefault_runMCMC = false;
		}else if (varName[i] == "useGlobalSamplingProbability"){
			_useGlobalSamplingProbability = stringToBool(varValue[i].c_str());
			isDefault_useGlobalSamplingProbability = false;
		}else if (varName[i] == "sampleProbsFilename"){
			_sampleProbsFilename = varValue[i];
			isDefault_sampleProbsFilename = false;
			//cout << left << setw(ppw) << "sampleProbsFilename" << "\t" << _sampleProbsFilename << endl;
	
		}else if (varName[i] == "globalSamplingFraction"){
			_globalSamplingFraction = atof(varValue[i].c_str());
			isDefault_globalSamplingFraction = false;
			//cout << left << setw(ppw) << "globalSamplingFraction" << "\t" << _globalSamplingFraction << endl;
		
		}else if (varName[i] == "updateLambdaInitScale"){
			_updateLambdaInitScale = atof(varValue[i].c_str());
			isDefault_updateLambdaInitScale = false;
			//cout << left << setw(ppw) << "updateLambdaInitScale" << "\t" << _updateLambdaInitScale << endl;
			
		}else if (varName[i] == "updateMuInitScale"){
			_updateMuInitScale = atof(varValue[i].c_str());
			isDefault_updateMuInitScale = false;
		}else if (varName[i] == "updateLambdaShiftScale"){
			_updateLambdaShiftScale = atof(varValue[i].c_str());
			isDefault_updateLambdaShiftScale = false;
		}else if (varName[i] == "updateMuShiftScale"){
			_updateMuShiftScale = atof(varValue[i].c_str());
			isDefault_updateMuShiftScale = false;
		}else if (varName[i] == "lambdaInit0"){
			_lambdaInit0 = atof(varValue[i].c_str());
			isDefault_lambdaInit0 = false;
		}else if (varName[i] == "lambdaShift0"){
			_lambdaShift0 = atof(varValue[i].c_str());
			isDefault_lambdaShift0 = false;
		}else if (varName[i] == "muInit0"){
			_muInit0 = atof(varValue[i].c_str());
			isDefault_muInit0 = false;
		}else if (varName[i] == "muShift0"){
			_muShift0 = atof(varValue[i].c_str());
			isDefault_muShift0 = false;
		}else if (varName[i] == "MeanSpeciationLengthFraction"){
			_MeanSpeciationLengthFraction = atof(varValue[i].c_str());
			isDefault_MeanSpeciationLengthFraction = false;
		}else if (varName[i] == "updateEventRateScale"){
			_updateEventRateScale = atof(varValue[i].c_str());
			isDefault_updateEventRateScale = false;
		}else if (varName[i] == "localGlobalMoveRatio"){
			_localGlobalMoveRatio = atof(varValue[i].c_str());	
			isDefault_localGlobalMoveRatio = false;
		}else if (varName[i] == "targetNumber"){
			_targetNumber = atof(varValue[i].c_str());
			isDefault_targetNumber = false;
		}else if (varName[i] == "lambdaInitPrior"){
			_lambdaInitPrior = atof(varValue[i].c_str());
			isDefault_lambdaInitPrior = false;
		}else if (varName[i] == "lambdaShiftPrior"){
			_lambdaShiftPrior= atof(varValue[i].c_str());
			isDefault_lambdaShiftPrior = false;
		}else if (varName[i] == "muInitPrior"){
			_muInitPrior = atof(varValue[i].c_str());
			isDefault_muInitPrior = false;
		}else if (varName[i] == "muShiftPrior"){
			_muShiftPrior = atof(varValue[i].c_str());
			isDefault_muShiftPrior = false;
		}else if (varName[i] == "segLength"){
			_segLength = atof(varValue[i].c_str());
			isDefault_segLength = false;
		}else if (varName[i] == "mcmcOutfile"){
			_mcmcOutfile = varValue[i];
			isDefault_mcmcOutfile = false;
		}else if (varName[i] == "eventDataOutfile"){
			_eventDataOutfile = varValue[i];
			isDefault_eventDataOutfile = false;
		}else if (varName[i] == "lambdaOutfile"){
			_lambdaOutfile = varValue[i];
			isDefault_lambdaOutfile = false;
		}else if (varName[i] == "muOutfile"){
			_muOutfile = varValue[i];
			isDefault_muOutfile = false;
		}else if (varName[i] == "acceptrateOutfile"){
			_acceptrateOutfile = varValue[i];
			isDefault_acceptrateOutfile = false;
		}else if (varName[i] == "lambdaNodeOutfile"){
			_lambdaNodeOutfile = varValue[i];
			isDefault_lambdaNodeOutfile = false;
		}else if (varName[i] == "treeWriteFreq"){
			_treeWriteFreq = atoi(varValue[i].c_str());
			isDefault_treeWriteFreq = false;
		}else if (varName[i] == "eventDataWriteFreq"){
			_eventDataWriteFreq = atoi(varValue[i].c_str());
			isDefault_eventDataWriteFreq = false;
		}else if (varName[i] == "mcmcWriteFreq"){
			_mcmcWriteFreq = atoi(varValue[i].c_str());
			isDefault_mcmcWriteFreq = false;
		}else if (varName[i] == "acceptWriteFreq"){
			_acceptWriteFreq = atoi(varValue[i].c_str());
			isDefault_acceptWriteFreq = false;
		}else if (varName[i] == "printFreq"){
			_printFreq = atoi(varValue[i].c_str());
			isDefault_printFreq = false;
		}else if (varName[i] == "NumberGenerations"){
			_NGENS = atoi(varValue[i].c_str());
			isDefault_NGENS = false;
		}else if (varName[i] == "updateRateEventNumber"){
			_updateRateEventNumber  = atof(varValue[i].c_str());
			isDefault_updateRateEventNumber = false;
		}else if (varName[i] == "updateRateEventPosition"){
			_updateRateEventPosition  = atof(varValue[i].c_str());
			isDefault_updateRateEventPosition = false;
		}else if (varName[i] == "updateRateEventRate"){
			_updateRateEventRate = atof(varValue[i].c_str());
			isDefault_updateRateEventRate = false;
		}else if (varName[i] == "updateRateLambda0"){
			_updateRateLambda0  = atof(varValue[i].c_str());
			isDefault_updateRateLambda0 = false;
		}else if (varName[i] == "updateRateLambdaShift"){
			_updateRateLambdaShift  = atof(varValue[i].c_str());
			isDefault_updateRateLambdaShift = false;
		}else if (varName[i] == "updateRateMu0"){
			_updateRateMu0  = atof(varValue[i].c_str());
			isDefault_updateRateMu0 = false;
		}else if (varName[i] == "updateRateMuShift"){
			_updateRateMuShift  = atof(varValue[i].c_str());
			isDefault_updateRateMuShift = false;
		}else if (varName[i] == "initialNumberEvents"){
			_initialNumberEvents  = atoi(varValue[i].c_str());
			isDefault_initialNumberEvents = false;
		}else if(varName[i] == "loadEventData" ){
			_loadEventData = stringToBool(varValue[i].c_str());	
			isDefault_loadEventData = false;
		}else if (varName[i] == "eventDataInfile"){
			_eventDataInfile = varValue[i].c_str();
			isDefault_eventDataInfile = false;
		}else if (varName[i] == "minCladeSizeForShift"){
			_minCladeSizeForShift = atoi(varValue[i].c_str());
			isDefault_minCladeSizeForShift = false;
		}else{
		// Parameter not found:
		//		add to list of potentially bad/misspelled params 
		//		and print for user.
			paramsNotFound.push_back(varName[i]);
		}	
		
	}

	cout << "Read a total of <<" << varName.size() << ">> parameter settings from control file" << endl;
	if (paramsNotFound.size() > 0){
		cout << endl << "********************************" << endl;
		cout << "BAMM error: one or more parameters from control file do not correspond" << endl;
		cout << "\tto valid model parameters.";
		cout << "Check the following to see if they are ";
		cout << "\tspecified (or spelled) correctly:" << endl << endl;
		for (int i = 0; i < paramsNotFound.size(); i++){
			 cout << setw(30) << paramsNotFound[i] << endl;		
		}
		cout << endl << "********************************" << endl << endl;
		cout << "Execution of BAMM terminated..." << endl;
		exit(1);
			 
	}
	
	//	Here we have a print block to output Settings:
	//	Any parameters NOT set will have the defaults.
	//	Thus user can specify a control file with ONLY 
	//		those parameters that they wish to change 
	//		from the defaults, eg inputfilename etc.

	//	Output list of default parameters.
	
	
}

bool Settings::stringToBool(const char * x){
	bool bval;
	if (atoi(x) == 0){
		bval = false;
	}else if (atoi(x) == 1){
		bval = true;
	}else{
		cout << "Invalid boolean when initializing class Settings" << endl;
		throw;
	}
	return bval;
}


void Settings::printCurrentSettings(bool printOnlyChangesToDefaults){

	int ppw = 29;
	

	
	cout << "*****************************************************" << endl;
	cout << "Current parameter settings: " << endl;
	if (printOnlyChangesToDefaults){
		cout << "\tPrinting only changes to default settings" << endl;
		cout << endl;
		if (!isDefault_treefile){
			cout << right << setw(ppw) << "treefile" << "\t\t" << _treefile << endl;
		}
		if (!isDefault_loadEventData){
			cout << right << setw(ppw) << "loadEventData" << "\t\t" << _loadEventData << endl;			
		}
		if (!isDefault_eventDataInfile){
			cout << right << setw(ppw) << "eventDataInfile" << "\t\t" << _eventDataInfile << endl;			
		}		
		if (!isDefault_sampleFromPriorOnly){
			cout << right << setw(ppw) << "sampleFromPriorOnly" << "\t\t" << _sampleFromPriorOnly << endl;
		}
		if (!isDefault_runMCMC){
			cout << right << setw(ppw) << "runMCMC" << "\t\t" << _runMCMC << endl;
		}
		if (!isDefault_initializeModel){
			cout << right << setw(ppw) << "initializeModel" << "\t\t" << _initializeModel << endl;
		}
		if (!isDefault_useGlobalSamplingProbability){
			cout << right << setw(ppw) << "useGlobalSamplingProbability" << "\t\t" << _useGlobalSamplingProbability << endl;
		}
		if (!isDefault_sampleProbsFilename){
			cout << right << setw(ppw) << "sampleProbsFilename" << "\t\t" << _sampleProbsFilename << endl;
		}
		if (!isDefault_globalSamplingFraction){
			cout << right << setw(ppw) << "globalSamplingFraction" << "\t\t" << _globalSamplingFraction << endl;
		}
		if (!isDefault_updateLambdaInitScale){
			cout << right << setw(ppw) << "updateLambdaInitScale" << "\t\t" << _updateLambdaInitScale << endl;
		}
		if (!isDefault_updateMuInitScale){
			cout << right << setw(ppw) << "updateMuInitScale" << "\t\t" << _updateMuInitScale << endl;
		}
		if (!isDefault_updateLambdaShiftScale){
			cout << right << setw(ppw) << "updateLambdaShiftScale" << "\t\t" << _updateLambdaShiftScale << endl;
		}
		if (!isDefault_updateMuShiftScale){
			cout << right << setw(ppw) << "updateMuShiftScale" << "\t\t" << _updateMuShiftScale << endl;
		}
		if (!isDefault_lambdaInit0){
			cout << right << setw(ppw) << "lambdaInit0" << "\t\t" << _lambdaInit0 << endl;
		}
		if (!isDefault_lambdaShift0){
			cout << right << setw(ppw) << "lambdaShift0" << "\t\t" << _lambdaShift0 << endl;
		}
		if (!isDefault_muInit0){
			cout << right << setw(ppw) << "muInit0" << "\t\t" << _muInit0 << endl;
		}
		if (!isDefault_muShift0){
			cout << right << setw(ppw) << "muShift0" << "\t\t" << _muShift0 << endl;
		}
		if (!isDefault_updateEventRateScale){
			cout << right << setw(ppw) << "updateEventRateScale" << "\t\t" << _updateEventRateScale << endl;
		}
		if (!isDefault_localGlobalMoveRatio){
			cout << right << setw(ppw) << "localGlobalMoveRatio" << "\t\t" << _localGlobalMoveRatio << endl;
		}
		if (!isDefault_targetNumber){
			cout << right << setw(ppw) << "targetNumber" << "\t\t" << _targetNumber << endl;
		}
		if (!isDefault_lambdaInitPrior){
			cout << right << setw(ppw) << "lambdaInitPrior" << "\t\t" << _lambdaInitPrior << endl;
		}
		if (!isDefault_lambdaShiftPrior){
			cout << right << setw(ppw) << "lambdaShiftPrior" << "\t\t" << _lambdaShiftPrior << endl;
		}
		if (!isDefault_muInitPrior){
			cout << right << setw(ppw) << "muInitPrior" << "\t\t" << _muInitPrior << endl;
		}
		if (!isDefault_muShiftPrior){
			cout << right << setw(ppw) << "muShiftPrior" << "\t\t" << _muShiftPrior << endl;
		}
		if (!isDefault_MeanSpeciationLengthFraction){
			cout << right << setw(ppw) << "MeanSpeciationLengthFraction" << "\t\t" << _MeanSpeciationLengthFraction << endl;
		}
		if (!isDefault_segLength){
			cout << right << setw(ppw) << "segLength" << "\t\t" << _segLength << endl;
		}
		if (!isDefault_mcmcOutfile){
			cout << right << setw(ppw) << "mcmcOutfile" << "\t\t" << _mcmcOutfile << endl;
		}
		if (!isDefault_lambdaOutfile){
			cout << right << setw(ppw) << "lambdaOutfile" << "\t\t" << _lambdaOutfile << endl;
		}
		if (!isDefault_muOutfile){
			cout << right << setw(ppw) << "muOutfile" << "\t\t" << _muOutfile << endl;
		}
		if (!isDefault_acceptrateOutfile){
			cout << right << setw(ppw) << "acceptrateOutfile" << "\t\t" << _acceptrateOutfile << endl;
		}
		if (!isDefault_lambdaNodeOutfile){
			cout << right << setw(ppw) << "lambdaNodeOutfile" << "\t\t" << _lambdaNodeOutfile << endl;
		}
		if (!isDefault_treeWriteFreq){
			cout << right << setw(ppw) << "treeWriteFreq" << "\t\t" << _treeWriteFreq << endl;
		}
		if (!isDefault_mcmcWriteFreq){
			cout << right << setw(ppw) << "mcmcWriteFreq" << "\t\t" << _mcmcWriteFreq << endl;
		}
		if (!isDefault_eventDataWriteFreq){
			cout << right << setw(ppw) << "eventDataWriteFreq" << "\t\t" << _eventDataWriteFreq << endl;
		}
		if (!isDefault_acceptWriteFreq){
			cout << right << setw(ppw) << "acceptWriteFreq" << "\t\t" << _acceptWriteFreq << endl;
		}
		if (!isDefault_printFreq){
			cout << right << setw(ppw) << "printFreq" << "\t\t" << _printFreq << endl;
		}
		if (!isDefault_NGENS){
			cout << right << setw(ppw) << "NGENS" << "\t\t" << _NGENS << endl;
		}
		if (!isDefault_updateRateEventNumber){
			cout << right << setw(ppw) << "updateRateEventNumber" << "\t\t" << _updateRateEventNumber << endl;
		}
		if (!isDefault_updateRateEventPosition){
			cout << right << setw(ppw) << "updateRateEventPosition" << "\t\t" << _updateRateEventPosition << endl;
		}
		if (!isDefault_updateRateEventRate){
			cout << right << setw(ppw) << "updateRateEventRate" << "\t\t" << _updateRateEventRate << endl;
		}
		if (!isDefault_updateRateLambda0){
			cout << right << setw(ppw) << "updateRateLambda0" << "\t\t" << _updateRateLambda0 << endl;
		}
		if (!isDefault_updateRateLambdaShift){
			cout << right << setw(ppw) << "updateRateLambdaShift" << "\t\t" << _updateRateLambdaShift << endl;
		}
		if (!isDefault_updateRateMu0){
			cout << right << setw(ppw) << "updateRateMu0" << "\t\t" << _updateRateMu0 << endl;
		}
		if (!isDefault_updateRateMuShift){
			cout << right << setw(ppw) << "updateRateMuShift" << "\t\t" << _updateRateMuShift << endl;
		}
		if (!isDefault_initialNumberEvents){
			cout << right << setw(ppw) << "initialNumberEvents" << "\t\t" << _initialNumberEvents << endl;
		}
		if (!isDefault_minCladeSizeForShift){
			cout << right << setw(ppw) << "minCladeSizeForShift" << "\t\t" << _minCladeSizeForShift << endl;
		}
		
	}else{
		cout << "\tPrinting ALL parameter settings " << endl;
		cout << endl;
		cout << right << setw(ppw) << "treefile" << "\t\t" << _treefile << endl;
		cout << right << setw(ppw) << "sampleFromPriorOnly" << "\t\t" << _sampleFromPriorOnly << endl;
		cout << right << setw(ppw) << "runMCMC" << "\t\t" << _runMCMC << endl;
		cout << right << setw(ppw) << "initializeModel" << "\t\t" << _initializeModel << endl;
		cout << right << setw(ppw) << "useGlobalSamplingProbability" << "\t\t" << _useGlobalSamplingProbability << endl;
		cout << right << setw(ppw) << "sampleProbsFilename" << "\t\t" << _sampleProbsFilename << endl;
		cout << right << setw(ppw) << "eventDataInfile" << "\t\t" << _eventDataInfile << endl;
		cout << right << setw(ppw) << "globalSamplingFraction" << "\t\t" << _globalSamplingFraction << endl;
		cout << right << setw(ppw) << "updateLambdaInitScale" << "\t\t" << _updateLambdaInitScale << endl;
		cout << right << setw(ppw) << "updateMuInitScale" << "\t\t" << _updateMuInitScale << endl;
		cout << right << setw(ppw) << "updateLambdaShiftScale" << "\t\t" << _updateLambdaShiftScale << endl;
		cout << right << setw(ppw) << "updateMuShiftScale" << "\t\t" << _updateMuShiftScale << endl;
		cout << right << setw(ppw) << "lambdaInit0" << "\t\t" << _lambdaInit0 << endl;
		cout << right << setw(ppw) << "lambdaShift0" << "\t\t" << _lambdaShift0 << endl;
		cout << right << setw(ppw) << "muInit0" << "\t\t" << _muInit0 << endl;
		cout << right << setw(ppw) << "muShift0" << "\t\t" << _muShift0 << endl;
		cout << right << setw(ppw) << "updateEventRateScale" << "\t\t" << _updateEventRateScale << endl;
		cout << right << setw(ppw) << "localGlobalMoveRatio" << "\t\t" << _localGlobalMoveRatio << endl;
		cout << right << setw(ppw) << "targetNumber" << "\t\t" << _targetNumber << endl;
		cout << right << setw(ppw) << "lambdaInitPrior" << "\t\t" << _lambdaInitPrior << endl;
		cout << right << setw(ppw) << "lambdaShiftPrior" << "\t\t" << _lambdaShiftPrior << endl;
		cout << right << setw(ppw) << "muInitPrior" << "\t\t" << _muInitPrior << endl;
		cout << right << setw(ppw) << "muShiftPrior" << "\t\t" << _muShiftPrior << endl;
		cout << right << setw(ppw) << "MeanSpeciationLengthFraction" << "\t\t" << _MeanSpeciationLengthFraction << endl;
		cout << right << setw(ppw) << "segLength" << "\t\t" << _segLength << endl;
		cout << right << setw(ppw) << "mcmcOutfile" << "\t\t" << _mcmcOutfile << endl;
		cout << right << setw(ppw) << "lambdaOutfile" << "\t\t" << _lambdaOutfile << endl;
		cout << right << setw(ppw) << "muOutfile" << "\t\t" << _muOutfile << endl;
		cout << right << setw(ppw) << "acceptrateOutfile" << "\t\t" << _acceptrateOutfile << endl;
		cout << right << setw(ppw) << "lambdaNodeOutfile" << "\t\t" << _lambdaNodeOutfile << endl;
		cout << right << setw(ppw) << "treeWriteFreq" << "\t\t" << _treeWriteFreq << endl;
		cout << right << setw(ppw) << "eventDataWriteFreq" << "\t\t" << _eventDataWriteFreq << endl;
		cout << right << setw(ppw) << "mcmcWriteFreq" << "\t\t" << _mcmcWriteFreq << endl;
		cout << right << setw(ppw) << "acceptWriteFreq" << "\t\t" << _acceptWriteFreq << endl;
		cout << right << setw(ppw) << "printFreq" << "\t\t" << _printFreq << endl;
		cout << right << setw(ppw) << "NGENS" << "\t\t" << _NGENS << endl;
		cout << right << setw(ppw) << "updateRateEventNumber" << "\t\t" << _updateRateEventNumber << endl;
		cout << right << setw(ppw) << "updateRateEventPosition" << "\t\t" << _updateRateEventPosition << endl;
		cout << right << setw(ppw) << "updateRateEventRate" << "\t\t" << _updateRateEventRate << endl;
		cout << right << setw(ppw) << "updateRateLambda0" << "\t\t" << _updateRateLambda0 << endl;
		cout << right << setw(ppw) << "updateRateLambdaShift" << "\t\t" << _updateRateLambdaShift << endl;
		cout << right << setw(ppw) << "updateRateMu0" << "\t\t" << _updateRateMu0 << endl;
		cout << right << setw(ppw) << "updateRateMuShift" << "\t\t" << _updateRateMuShift << endl;
		cout << right << setw(ppw) << "initialNumberEvents" << "\t\t" << _initialNumberEvents << endl;
		cout << right << setw(ppw) << "minCladeSizeForShift" << "\t\t" << _minCladeSizeForShift << endl;
		
	}
	
	
	cout << endl;
	
	cout << "*****************************************************" << endl;	
	
	
	
	
}


/*
 
 This (or another function) needs to be re-written to deal with potential input files 
 for TRAIT analysis versus SPECIATION-EXTINCTION analysis.
 
 */

void Settings::parseCommandLineInput(int argc, vector<string> &instrings, string modeltype){
	
	
	vector<string> badFlags;
	int i = 0;
	while (i + 1 != argc){
		if (instrings[i] == "-control"){
			string controlfile = instrings[i + 1];
			cout << "\n\nInitializing BAMM using control file <<" << controlfile << ">>" << endl;
			ifstream instream(controlfile.c_str());
			if (!instream){
				cout << "File not found error: cannot locate control file in\n";
				cout <<	 "specified directory. Exiting BAMM" << endl << endl;
				exit(1);
			}else{
				if (modeltype == "speciationextinction"){
					initializeSettings(controlfile);				
				}else if (modeltype == "trait"){
					trait_initializeSettings(controlfile);
				}else{
					cout << "Unsupported modeltype in Settings::parseCommandLineInput" << endl;
				}
		
				break; // exit for loop if initialized with controlfile...
			}
		}else{
			if (i != 0){
				// If i == 0, we just assume that this is either (i) the path to the executable
				//	or (ii), the commands used to invoke it.
				//	May depend on system/compiler, but first argument should
				//	not be a valid flag...
				badFlags.push_back(instrings[i]);
			}
			i++;
		}		
	
	}
	
	if (areAllParametersSetToDefaults()){
		cout << "Failed to initialize parameter values\nExiting BAMM" << endl;
		exit(1);
	
	}

}


 




