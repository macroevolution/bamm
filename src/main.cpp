
/*

 This version of BAMM does speciation-extinction and trait evolution.

*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Tree.h"
#include "MbRandom.h"
#include "Model.h"
#include "MCMC.h"
#include "Settings.h"
#include "TraitMCMC.h"
#include "TraitModel.h"
#include "Autotune.h"

void usage () {
    std::cout << std::endl << "Program usage:" << std::endl;
    std::cout << "./bamm -control control_filename" << std::endl<< std::endl;
}
  
int main (int argc, char* argv[])
{

    if (argc == 1) {
        usage();
        exit(0);
    } else if (argc == 2) {
        std::string arg = std::string(argv[1]);
        if (arg == "-h" || arg == "h" || arg == "help" || arg == "-help") {
            usage();
            exit(0);
        } else {
            std::cout << "Unknown argument '" << arg << "'." << std::endl;
            usage();
            exit(0);
        }
    }
    
//    MbRandom myRNG; // *** don't do this until AFTER reading control file ***
    Settings mySettings;
    for (int i = 0; i < 20; i++) {
        std::cout << "#";
    }
    
    if (argc <= 1) {
        std::cout << "Removed option to initialized BAMM with default settings" << std::endl;
        std::cout << "You must specify a controlfilename" << std::endl;
        throw;
        
    } else if (argc > 1) {
        // IF > 1 things read assume other args:

        std::vector<std::string> instrings;
        for (int i = 0; i < argc; i++) {
            instrings.push_back(argv[i]);
        }
        
        mySettings.parseCommandLineInput(argc, instrings);
    
    } else {
        std::cout << "Uninterpretable input. Exiting BAMM." << std::endl;
        exit(1);
    }
    
    MbRandom myRNG(mySettings.getSeed());
    std::cout << mySettings.getModeltype() << std::endl;
    
/*
    double terp = 0.0;
    std::cout << "TESTING uniformRv:" << std::endl;
    for (int i = 0; i < 10000000; i++) {
        terp = myRNG.uniformRv();
        if (i > 9999975) { // only print out the last few numbers
            std::cout << i << ": " << terp << std::endl;
        }
    }
    exit(0);
*/
    
    if (mySettings.getModeltype() == "speciationextinction") {
        std::cout << "Initializing diversification (speciationextinction) model" << std::endl;
        for (int i = 0; i < 20; i++) {
            std::cout << "#";
        }
        
        std::cout << std::endl << std::endl  << "SPECIATION-EXTINCTION BAMM" << std::endl << std::endl;
    
        mySettings.printCurrentSettings();
        mySettings.checkSettingsAreUserDefined();
        std::string treefile = mySettings.getTreeFilename();
        Tree intree(treefile, &myRNG);
        
        if (mySettings.getUseGlobalSamplingProbability()) {
            std::cout << "Initializing with global sampling probability\n" << std::endl;
            intree.initializeSpeciationExtinctionModel(mySettings.getGlobalSamplingFraction());
        } else {
            std::cout << "Species-specific sampling fractions are not validated yet...\n" <<
            std::endl;
            // code should be supported for this but need to check..
            intree.initializeSpeciationExtinctionModel(mySettings.getSampleProbsFilename());
            //throw;
        }
        
        intree.setCanNodeHoldEventByDescCount(mySettings.getMinCladeSizeForShift());
        intree.setTreeMap(intree.getRoot());
        if (mySettings.getInitializeModel() && !mySettings.getRunMCMC()) {
            Model myModel(&myRNG, &intree, &mySettings);
            std::cout << "Initializing model but not running MCMC" << std::endl;
        } else if (mySettings.getInitializeModel() && mySettings.getAutotune()) {
            Model myModel(&myRNG, &intree, &mySettings);
            Autotune myTuneObject(&myRNG, &myModel, &mySettings);
        } else if (mySettings.getInitializeModel() && mySettings.getRunMCMC()) {
            Model myModel(&myRNG, &intree, &mySettings);
            MCMC myMCMC(&myRNG, &myModel, &mySettings);
        } else {
            std::cout << "Unsupported option in main....\n" << std::endl;
        }
        
    } else if (mySettings.getModeltype() == "trait") {
        std::cout << "Initializing phenotypic evolution (trait) model " << std::endl;
        for (int i = 0; i < 20; i++) {
            std::cout << "#";
        }
        std::cout << std::endl << std::endl  << "TRAIT DIVERSIFICATION BAMM" << std::endl << std::endl;
        
        mySettings.printCurrentSettings();
        mySettings.checkSettingsAreUserDefined();

        std::string treefile = mySettings.getTreeFilename();
        Tree intree(treefile, &myRNG);
        
        intree.setAllNodesCanHoldEvent();
        intree.setTreeMap(intree.getRoot());

        intree.getPhenotypesMissingLatent(mySettings.getTraitFile());
        intree.initializeTraitValues();
        
        if (mySettings.getInitializeModel() && !mySettings.getRunMCMC()) {
            std::cout << "Initializing model but not running MCMC" << std::endl;
            TraitModel myModel(&myRNG, &intree, &mySettings);
        } else if (mySettings.getInitializeModel() && mySettings.getRunMCMC() && !mySettings.getAutotune()) {
            std::cout << "Initializing model and MCMC chain" << std::endl;
            TraitModel myModel(&myRNG, &intree, &mySettings);
            TraitMCMC myMCMC(&myRNG, &myModel, &mySettings);
        } else if (mySettings.getInitializeModel() && mySettings.getAutotune()) {
            std::cout << "Autotune option not yet supported for phenotypic (trait) analysis" << std::endl;
            exit(0);
        } else {
            std::cout << "Invalid run settings specified in main\n" << std::endl;
            exit(0);
        }

    } else if (mySettings.getModeltype() == "EMPTY_STRING") {
        std::cout << "You did not specify a modeltype." << std::endl;
        std::cout << "You must specify one of the following for parameter modeltype\n\n" << std::endl;
        std::cout << "\t\tspeciationextinction\n\t\ttrait\n\n" << std::endl << std::endl;
    } else {
        std::cout << "Invalid modeltype specification in controlfile" << std::endl;
        std::cout << "Specified modeltype was <<" << mySettings.getModeltype() << ">> \n" << std::endl;
        std::cout << "Valid options: \n" << "\t\tspeciationextinction\n\t\ttrait\n\n" << std::endl;
    }

    return 0;
}
