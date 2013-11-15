
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


int main (int argc, char* argv[])
{

    //argc =3;
    //argv[1] = "-control";
    //argv[2] = "control.txt";


    //for (int i = 0; i < argc; i++){
    //  std::cout << argc << "\t" << argv[i] << std::endl;
    //}


    std::string modeltype = std::string(argv[1]);
    //std::string modeltype = "trait";

    MbRandom myRNG;
    Settings mySettings;

    if (modeltype == "speciationextinction") {

        for (int i = 0; i < 20; i++)
            std::cout << "#";

        std::cout << std::endl << std::endl  << "SPECIATION-EXTINCTION BAMM" << std::endl << std::endl;

        for (int i = 0; i < 20; i++)
            std::cout << "#";
        if (argc <= 1) {
            std::cout << "\nInitializing BAMM with default settings." << std::endl;
            std::cout << "\tThis may not be OK - consult manual for usage information\n" << std::endl;
            mySettings.initializeSettings();
        } else if (argc > 1) {
            // IF > 1 things read assume other args:
            std::vector<std::string> instrings;
            for (int i = 0; i < argc; i++)
                if (i != 1) // skip trait or speciationextinction option
                    instrings.push_back(argv[i]);
            argc--;
            mySettings.parseCommandLineInput(argc, instrings, modeltype);
        } else {
            std::cout << "Uninterpretable input. Exiting BAMM." << std::endl;
            exit(1);
        }


        mySettings.printCurrentSettings(true);
        std::string treefile = mySettings.getTreeFilename();
        Tree intree(treefile, &myRNG);

        if (mySettings.getUseGlobalSamplingProbability()) {
            std::cout << "Initializing with global sampling probability\n" << std::endl;
            intree.initializeSpeciationExtinctionModel(
                mySettings.getGlobalSamplingFraction());

        } else {
            std::cout << "Species-specific sampling fractions are not validated yet...\n" <<
                 std::endl;
            // code should be supported for this but need to check..
            intree.initializeSpeciationExtinctionModel(mySettings.getSampleProbsFilename());
            //throw;
        }

        //intree.printCanHoldEventByNode();


        std::cout << std::endl << std::endl;
        //intree.setAllNodesCanHoldEvent();

        std::cout << "MinCladeSize: " << mySettings.getMinCladeSizeForShift() << std::endl;

        intree.setCanNodeHoldEventByDescCount(mySettings.getMinCladeSizeForShift());

        intree.setTreeMap(intree.getRoot());

        //intree.printCanHoldEventByNode();

        if (mySettings.getInitializeModel() && !mySettings.getRunMCMC()) {
            Model myModel(&myRNG, &intree, &mySettings);
            std::cout << "Initializing model but not running MCMC" << std::endl;

        } else if (mySettings.getInitializeModel() && mySettings.getRunMCMC()) {
            Model myModel(&myRNG, &intree, &mySettings);
            MCMC myMCMC(&myRNG, &myModel, &mySettings);

        } else
            std::cout << "Unsupported option in main....\n" << std::endl;





    } else if (modeltype == "trait") {

        for (int i = 0; i < 20; i++)
            std::cout << "#";

        std::cout << std::endl << std::endl  << "TRAIT BAMM" << std::endl << std::endl;

        for (int i = 0; i < 20; i++)
            std::cout << "#";

        if (argc <= 1) {
            std::cout << "\nInitializing BAMMt with default settings." << std::endl;
            std::cout << "\tThis may not be OK - consult manual for usage information\n" << std::endl;
            mySettings.trait_initializeSettings();
        } else if (argc > 1) {
            // IF > 1 things read assume other args:
            std::vector<std::string> instrings;
            for (int i = 0; i < argc; i++)
                if (i != 1) // skip trait or speciationextinction option
                    instrings.push_back(argv[i]);
            argc--;


            mySettings.parseCommandLineInput(argc, instrings, modeltype);
            mySettings.checkAreTraitInitialSettingsValid();
        } else {
            std::cout << "Uninterpretable input. Exiting BAMM." << std::endl;
            exit(1);
        }


        //mySettings.trait_printCurrentSettings(true);
        std::string treefile = mySettings.getTreeFilename();
        Tree intree(treefile, &myRNG);

        intree.setAllNodesCanHoldEvent();
        intree.setTreeMap(intree.getRoot());

        //intree.getPhenotypes(mySettings.getTraitFile());
        intree.getPhenotypesMissingLatent(mySettings.getTraitFile());


        intree.initializeTraitValues();


        if (mySettings.getInitializeModel() && !mySettings.getRunMCMC()) {
            std::cout << "Initializing model but not running MCMC" << std::endl;
            TraitModel myModel(&myRNG, &intree, &mySettings);
        }

        if (mySettings.getInitializeModel() && mySettings.getRunMCMC()) {
            std::cout << "Initializing model and MCMC chain" << std::endl;
            TraitModel myModel(&myRNG, &intree, &mySettings);
            TraitMCMC myMCMC(&myRNG, &myModel, &mySettings);

            //intree.echoMeanBranchTraitRates();

        }


    } else {
        std::cout << "Unsupported analysis" << std::endl;
        exit(1);
    }




    return 0;
}




