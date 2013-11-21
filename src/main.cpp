#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <cstdlib>

#include "Tree.h"
#include "MbRandom.h"
#include "Model.h"
#include "MCMC.h"
#include "Settings.h"
#include "TraitMCMC.h"
#include "TraitModel.h"
#include "Autotune.h"


const char* currentTime();
void exitWithMessageUsage();
void exitWithErrorUnknownArgument(const std::string& arg);
void exitWithErrorNoControlFile();

  
int main (int argc, char* argv[])
{
    if (argc == 1) {
        exitWithMessageUsage();
    }

    std::string controlFilename;

    // Process command-line arguments
    int i = 1;
    while (i < argc) {
        std::string arg = std::string(argv[i]);
        if (arg == "-h" || arg == "--help" || arg == "-help") {
            exitWithMessageUsage();
        } else if (arg == "-c" || arg == "--control" || arg == "-control") {
            if (++i < argc) {
                controlFilename = std::string(argv[i]);
            } else {
                exitWithErrorNoControlFile();
            }
        } else {
            exitWithErrorUnknownArgument(arg);
        }

        i++;
    }

    // Load settings from control file
    Settings mySettings(controlFilename);

    for (int i = 0; i < 20; i++) {
        std::cout << "#";
    }
    std::cout << "\n";
    
    MbRandom myRNG(mySettings.getSeed());

    std::cout << "Random seed:\t\t" << myRNG.getSeed() << std::endl;

    std::string commandLine(argv[0]);
    for (int i = 1; i < argc; i++) {
        commandLine += std::string(" ") + argv[i];
    }

    std::ofstream runInfoFile(mySettings.getRunInfoFilename().c_str());
    runInfoFile << "command line: " << commandLine << "\n";
    runInfoFile << "git commit id: " << GIT_COMMIT_ID << "\n";
    runInfoFile << "random seed: " << myRNG.getSeed() << "\n";
    runInfoFile << "start time: " << currentTime() << "\n";

    std::cout << mySettings.getModeltype() << std::endl;
    
    if (mySettings.getModeltype() == "speciationextinction") {
        std::cout << "Initializing diversification (speciationextinction) " <<
            "model\n";

        for (int i = 0; i < 20; i++) {
            std::cout << "#";
        }
        
        std::cout << "\nSPECIATION-EXTINCTION BAMM\n\n";
    
        mySettings.printCurrentSettings();
        mySettings.printCurrentSettings(runInfoFile);

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
        mySettings.printCurrentSettings(runInfoFile);

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
    }

    runInfoFile << "end time: " << currentTime();
    runInfoFile.close();

    return 0;
}

const char* currentTime()
{
    time_t curTime;
    time(&curTime);
    return std::ctime(&curTime);
}


void exitWithMessageUsage()
{
    std::cout << "Usage: ./bamm -c control_filename\n";
    std::exit(0);
}


void exitWithErrorUnknownArgument(const std::string& arg)
{
    std::cout << "Unknown argument " << arg << ".\n";
    std::exit(1);
}


void exitWithErrorNoControlFile()
{
    std::cout << "ERROR: No control file specified.\n";
    std::cout << "Fix by specifying a control file name.\n";
    std::exit(1);
}
