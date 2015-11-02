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


void printAboutInformation();
std::string buildCommandLine(int argc, char* argv[]);
const char* currentTime();
ModelFactory* createModelFactory(const std::string& modelType);


int main (int argc, char* argv[])
{
    // Process command-line arguments and load settings
    CommandLineProcessor commandLine(argc, argv);
    Settings settings(commandLine.controlFileName(), commandLine.parameters());

    printAboutInformation();

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
        
        std::cout << "\n***********************************************" << std::endl;
        std::cout << "You selected to simulate the prior distribution of shifts." << std::endl;
        std::cout << "This option has been disabled, as the " << std::endl;
        std::cout << "exact (analytical) prior is now implemented in BAMMtools. " << std::endl;
        std::cout << "Please consult the website or BAMMtools documentation " << std::endl;
        std::cout << "for more information\n\n" << std::endl;
        //FastSimulatePrior fsp(random, &settings);
    }

    log(Message, runInfoFile) << "End time: " << currentTime() << "\n";
    runInfoFile.close();

    return 0;
}


void printAboutInformation()
{
    log(Message) << "BAMM " << BAMM_VERSION
        << " (" << BAMM_VERSION_DATE << ")\n"
        << "Copyright (C) 2012-2014 Daniel Rabosky\n"
        << "BAMM is distributed under the GNU General Public License.\n"
        << "See http://bamm-project.org for more information.\n\n";
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
