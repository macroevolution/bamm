#include "CommandLineProcessor.h"
#include "Settings.h"
#include "Random.h"
#include "ModelFactory.h"
#include "SpExModelFactory.h"
#include "TraitModelFactory.h"
#include "FastSimulatePrior.h"
#include "MetropolisCoupledMCMC.h"
#include "Log.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <cstdlib>


void printAboutBox();
std::string buildCommandLine(int argc, char* argv[]);
const char* currentTime();
ModelFactory* createModelFactory(const std::string& modelType);


int main (int argc, char* argv[])
{
    printAboutBox();

    // Process command-line arguments and load settings
    CommandLineProcessor commandLine(argc, argv);
    Settings settings(commandLine.controlFileName(), commandLine.parameters());

    // Setup pseudorandom generator
    int seed = settings.get<long int>("seed");
    Random random = (seed > 0) ? Random(seed) : Random();
    seed = random.getSeed();    // Get actual seed in case it is based on clock
    log(Message) << "Random seed: " << seed << "\n";

    // Setup "run info" file and print current settings
    std::ofstream runInfoFile(settings.get("runInfoFilename").c_str());
    log(Message, runInfoFile) << "Command line: "
        << buildCommandLine(argc, argv) << "\n";
    log(Message, runInfoFile) << "Git commit id: " << GIT_COMMIT_ID << "\n";
    log(Message, runInfoFile) << "Random seed: " << seed << "\n";
    log(Message, runInfoFile) << "Start time: " << currentTime() << "\n";
    settings.printCurrentSettings(runInfoFile);

    // Create model factory based on model type
    ModelFactory* modelFactory = createModelFactory(settings.get("modeltype"));

    if (settings.get<bool>("initializeModel")) {
        // MetropolisCoupledMCMC will initialize model(s)
        MetropolisCoupledMCMC mc3(random, settings, modelFactory);

        if (settings.get<bool>("runMCMC")) {
            mc3.run();
        }
    }

    delete modelFactory;

    if (settings.get<bool>("simulatePriorShifts")){
        FastSimulatePrior fsp(random, &settings);
    }

    log(Message, runInfoFile) << "End time: " << currentTime() << "\n";
    runInfoFile.close();

    return 0;
}


void printAboutBox()
{
    log(Message) << "\
+----------------------------------------------------------------------+\n\
|   BAMM: Bayesian Analysis of Macroevolutionary Mixtures              |\n\
+----------------------------------------------------------------------+\n\
|                                                                      |\n\
|   Daniel Rabosky <drabosky@umich.edu>                                |\n\
|   University of Michigan, Ann Arbor, MI, USA                         |\n\
|                                                                      |\n\
|   Authors: Carlos Anderson, Joseph Brown, Michael Grundler,          |\n\
|            Daniel Rabosky, Jeff Shi, Pascal Title                    |\n\
|                                                                      |\n\
|   Copyright (C) 2012-2014 Daniel Rabosky                             |\n\
|   See LICENSE for details.                                           |\n\
|                                                                      |\n\
+----------------------------------------------------------------------+\n\n";
}


std::string buildCommandLine(int argc, char* argv[])
{
    // String together the command-line arguments
    std::string commandLine(argv[0]);
    for (int i = 1; i < argc; i++) {
        commandLine += std::string(" ") + argv[i];
    }

    return commandLine;
}


const char* currentTime()
{
    time_t curTime;
    time(&curTime);
    return std::ctime(&curTime);
}


ModelFactory* createModelFactory(const std::string& modelType)
{
    if (modelType == "speciationextinction") {
        log(Message) << "\nModel type: Speciation/Extinction\n";
        return new SpExModelFactory();
    } else if (modelType == "trait") {
        log(Message) << "\nModel type: Trait\n";
        return new TraitModelFactory();
    } else {
        exitWithError("Unrecognized model");
        return NULL;    // Never reached but suppresses warning
    }
}
