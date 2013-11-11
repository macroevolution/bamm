
/* 
 
 This version of BAMM does speciation-extinction and trait evolution.
 
 
 */



#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "event.h"
#include "node.h"
#include "MbRandom.h"
#include "model.h"
#include "BranchHistory.h"
#include "mcmc.h"
#include "Settings.h"
#include "TraitMCMC.h"
#include "TraitModel.h"


// Specify the type of BAMM run at compile-time
//#define BUILD_TYPE_DEFAULT "speciationextinction"
#define BUILD_TYPE_DEFAULT "trait"



int main (int argc, char * argv[]) {

	//argc =3;
	//argv[1] = "-control";
	//argv[2] = "control.txt";
	
 
	//for (int i = 0; i < argc; i++){
	//	cout << argc << "\t" << argv[i] << endl;
	//}

	
	string modeltype = BUILD_TYPE_DEFAULT;
 	//string modeltype = "trait";
 	
 	MbRandom myRNG;
	Settings mySettings;
	
	if (modeltype == "speciationextinction"){
		
		for (int i = 0; i < 20; i++){
			cout << "#";
		}
		
		cout << endl << endl  << "SPECIATION-EXTINCTION BAMM" << endl << endl;
		
		for (int i = 0; i < 20; i++){
			cout << "#";
		}		
		if (argc <= 1){
			cout << "\nInitializing BAMM with default settings." << endl;
			cout << "\tThis may not be OK - consult manual for usage information\n" << endl;
			mySettings.initializeSettings();
		}else if (argc > 1){
			// IF > 1 things read assume other args:
			vector<string> instrings;
			for (int i = 0; i < argc; i++){
				instrings.push_back(argv[i]);
			}
			mySettings.parseCommandLineInput(argc, instrings, modeltype);	
		}else{
			cout << "Uninterpretable input. Exiting BAMM." << endl;
			exit(1);
		}
		
		
		mySettings.printCurrentSettings(true);
		string treefile = mySettings.getTreeFilename();
		Tree intree(treefile, &myRNG);	
		
		if (mySettings.getUseGlobalSamplingProbability()){
			cout << "Initializing with global sampling probability\n" << endl;
			intree.initializeSpeciationExtinctionModel(mySettings.getGlobalSamplingFraction());	
			
		}else{
			cout << "Species-specific sampling fractions are not validated yet...\n" << endl;
			// code should be supported for this but need to check..
			intree.initializeSpeciationExtinctionModel(mySettings.getSampleProbsFilename());		
			//throw;
		}
		
		//intree.printCanHoldEventByNode();		
 		
 		
 		cout << endl << endl;
		//intree.setAllNodesCanHoldEvent();
		
		cout << "MinCladeSize: " << mySettings.getMinCladeSizeForShift() << endl;
				
		intree.setCanNodeHoldEventByDescCount(mySettings.getMinCladeSizeForShift());
		
		intree.setTreeMap(intree.getRoot());

		//intree.printCanHoldEventByNode();

		if (mySettings.getInitializeModel() && !mySettings.getRunMCMC()){
			Model myModel(&myRNG, &intree, &mySettings);		
			cout << "Initializing model but not running MCMC" << endl;
			
		}else if (mySettings.getInitializeModel() && mySettings.getRunMCMC()){
			Model myModel(&myRNG, &intree, &mySettings);		
			MCMC myMCMC(&myRNG, &myModel, &mySettings);
		
		}else{
			cout << "Unsupported option in main....\n" << endl;
		}
		
		
		
		
		
	}else if (modeltype == "trait"){
		
		for (int i = 0; i < 20; i++){
			cout << "#";
		}
	
		cout << endl << endl  << "TRAIT BAMM" << endl << endl;

		for (int i = 0; i < 20; i++){
			cout << "#";
		}		
		
		if (argc <= 1){
			cout << "\nInitializing BAMMt with default settings." << endl;
			cout << "\tThis may not be OK - consult manual for usage information\n" << endl;
			mySettings.trait_initializeSettings();
		}else if (argc > 1){
			// IF > 1 things read assume other args:
			vector<string> instrings;
			for (int i = 0; i < argc; i++){
				instrings.push_back(argv[i]);
			}
			
			
			mySettings.parseCommandLineInput(argc, instrings, modeltype);
			mySettings.checkAreTraitInitialSettingsValid();
		}else{
			cout << "Uninterpretable input. Exiting BAMM." << endl;
			exit(1);
		}
		
		
		//mySettings.trait_printCurrentSettings(true);
		string treefile = mySettings.getTreeFilename();
		Tree intree(treefile, &myRNG);		
		
		intree.setAllNodesCanHoldEvent();
		intree.setTreeMap(intree.getRoot());
		
		//intree.getPhenotypes(mySettings.getTraitFile());
		intree.getPhenotypesMissingLatent(mySettings.getTraitFile());
		
		
		intree.initializeTraitValues();
		
		
		if (mySettings.getInitializeModel() && !mySettings.getRunMCMC()){
			cout << "Initializing model but not running MCMC" << endl;
			TraitModel myModel(&myRNG, &intree, &mySettings);				
		}
		
		if (mySettings.getInitializeModel() && mySettings.getRunMCMC()){
			cout << "Initializing model and MCMC chain" << endl;
			TraitModel myModel(&myRNG, &intree, &mySettings);		
			TraitMCMC myMCMC(&myRNG, &myModel, &mySettings);
		
			//intree.echoMeanBranchTraitRates();
		
		}
 

	}else {
		cout << "Unsupported analysis" << endl;
		exit(1);
	}	
	
	
	
	
    return 0;
}




