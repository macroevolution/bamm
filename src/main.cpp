#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <cstdlib>

#include "Tree.h"
#include "Node.h"
#include "MbRandom.h"
#include "MCMC.h"
#include "SpExModel.h"
#include "Settings.h"
#include "TraitModel.h"
#include "FastSimulatePrior.h"
#include "Prior.h"
#include "SpExDataWriter.h"
#include "TraitDataWriter.h"


void printAbout();
std::string buildCommandLine(int argc, char* argv[]);
const char* currentTime();
void exitWithMessageUsage();
void exitWithErrorUnknownArgument(const std::string& arg);
void exitWithErrorNoControlFile();

#include "Log.h"

int main (int argc, char* argv[])
{
    printAbout();

    if (argc == 1) {
        exitWithMessageUsage();
    }

    std::vector<UserParameter> commandLineParameters;
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
        // Let the Settings class (below) process the remaining arguments
        } else {
            if (++i < argc) {
                // First, make sure argument starts with "--"
                if ((arg[0] != '-') || (arg[1] != '-')) {
                    log(Error) << "Invalid command-line option <<"
                        << arg << ">>.\n";
                    std::exit(1);
                }
                arg = arg.substr(2);    // Cut out the starting "--"
                std::string arg_value = std::string(argv[i]);
                UserParameter param(arg, arg_value);
                commandLineParameters.push_back(param);
            } else {
                log(Error) << "Missing a command-line argument.\n";
                std::exit(1);
            }
        }

        i++;
    }

    // Load settings from control file
    Settings settings(controlFilename, commandLineParameters);

    int seed = settings.get<int>("seed");
    MbRandom rng(seed);
    seed = rng.getSeed();    // Get actual seed in case it is based on clock

    Prior prior(&rng, &settings);

    std::ofstream runInfoFile(settings.get("runInfoFilename").c_str());
    log(Message, runInfoFile) << "Command line: "
        << buildCommandLine(argc, argv) << "\n";
    //log(Message, runInfoFile) << "Git commit id: " << GIT_COMMIT_ID << "\n";
    log(Message, runInfoFile) << "Random seed: " << seed << "\n";
    log(Message, runInfoFile) << "Start time: " << currentTime() << "\n";
    
    log(Message) << "Random seed: " << seed << "\n";

    if (settings.get("modeltype") == "speciationextinction") {
        log(Message) << "\nModel type: Speciation/Extinction\n";

        settings.printCurrentSettings(runInfoFile);

        if (settings.get<bool>("initializeModel")) {
            SpExModel model(&rng, &settings, &prior);

            if (settings.get<bool>("runMCMC")) {
                int numberOfGenerations =
                    settings.get<int>("numberGenerations");
                SpExDataWriter dataWriter(settings, model);
                MCMC mcmc(rng, model, numberOfGenerations, dataWriter);
                mcmc.run();
            }
        }

    } else if (settings.get("modeltype") == "trait") {
        log(Message) << "\nModel type: Trait\n";
        
        settings.printCurrentSettings(runInfoFile);

        if (settings.get<bool>("initializeModel")) {
            TraitModel model(&rng, &settings, &prior);

            if (settings.get<bool>("runMCMC")) {
                int numberOfGenerations =
                    settings.get<int>("numberGenerations");
                TraitDataWriter dataWriter(settings, model);
                MCMC mcmc(rng, model, numberOfGenerations, dataWriter);
                mcmc.run();
            }
        }
    }

    if (settings.get<bool>("simulatePriorShifts")){
        FastSimulatePrior fsp(&rng, &settings);
    }

    log(Message, runInfoFile) << "End time: " << currentTime() << "\n";
    runInfoFile.close();

    return 0;
}


void printAbout()
{
    log(Message) << "\
+--------------------------------------------------------------------------+\n\
|   BAMM: Bayesian Analysis of Macroevolutionary Mixtures                  |\n\
+--------------------------------------------------------------------------+\n\
|                                                                          |\n\
|   Daniel Rabosky <drabosky@umich.edu>                                    |\n\
|   University of Michigan, Ann Arbor, MI, USA                             |\n\
|                                                                          |\n\
|   Authors: Carlos Anderson, Joseph Brown, Michael Grundler,              |\n\
|            Daniel Rabosky, Jeff Shi, Pascal Title                        |\n\
|                                                                          |\n\
|   Copyright (C) 2012-2014 Daniel Rabosky                                 |\n\
|   See LICENSE for details.                                               |\n\
|                                                                          |\n\
+--------------------------------------------------------------------------+\n\
\n";
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
