#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "MCMC.h"
#include "Model.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "Settings.h"
#include "BranchEvent.h"


MCMC::MCMC(MbRandom* ran, Model* mymodel, Settings* sp)
{
    std::cout << "Initializing MCMC object..." << std::endl;

    ranPtr = ran;
    ModelPtr = mymodel;
    sttings = sp;

    //check if file exists; delete;

    // MCMC parameters - for now....//
    _mcmcOutFilename       = sttings->getMCMCoutfile();
    _lambdaOutFilename     = sttings->getLambdaOutfile();
    _muOutFilename         = sttings->getMuOutfile();
    _acceptOutFilename     = sttings->getAcceptrateOutfile();
    _lambdaNodeOutFilename = sttings->getLambdaNodeOutfile();
    _eventDataOutFilename  = sttings->getEventDataOutfile();

    _treeWriteFreq =        sttings->getBranchRatesWriteFreq();
    _eventDataWriteFreq =   sttings->getEventDataWriteFreq();
    _mcmcWriteFreq =        sttings->getMCMCwriteFreq();
    _acceptWriteFreq =      sttings->getAcceptWriteFreq();
    _printFreq =            sttings->getPrintFreq();
    _NGENS =                sttings->getNGENS();

    _firstLine = true; // Print header for the first line of output
    
    bool fileOverwrite =     sttings->getOverwrite();

    if (!fileOverwrite) {
        if (anyOutputFileExists()) {
            exitWithErrorOutputFileExists();
        }
    }

    // Open streams for writing (overwrite)
    _mcmcOutStream.open(_mcmcOutFilename.c_str());
    _lambdaOutStream.open(_lambdaOutFilename.c_str());
    _muOutStream.open(_muOutFilename.c_str());
    _acceptOutStream.open(_acceptOutFilename.c_str());
    _lambdaNodeOutStream.open(_lambdaNodeOutFilename.c_str());
    _eventDataOutStream.open(_eventDataOutFilename.c_str());

    setUpdateWeights();
    ModelPtr->resetGeneration();

    for (int i = 0; i < sttings->getInitialNumberEvents(); i++) {
        ModelPtr->addEventToTree();
    }

    std::cout << "MCMC object successfully initialized." << std::endl << std::endl;
    std::cout << "Running MCMC chain for " << _NGENS << " generations." << std::endl << std::endl;

    std::cout << std::setw(10) << "Generation" << std::setw(10) << "lnLik" <<  std::setw(10);
    std::cout << "N_shifts" << std::setw(15) << "LogPrior" << std::setw(15) << "acceptRate" <<
         std::endl;

    /*run chain*/
    for (int i = 0; i < _NGENS; i++) {
        int parmToUpdate = pickParameterClassToUpdate();
        updateState(parmToUpdate); // update state

        if ((i % _treeWriteFreq) == 0) {

            writeBranchSpeciationRatesToFile();
            writeBranchExtinctionRatesToFile();
            // Deprecating this: no need to write this
            //writeNodeSpeciationRatesToFile();
        }

        if ((i % _eventDataWriteFreq) == 0) {
            writeEventDataToFile();
        }
        
        if ((i % _mcmcWriteFreq) == 0) {
            writeStateToFile();
        }

        if ((i % _printFreq == 0)) {
            printStateData();
        
            // Resets acceptance rates every time
            //  state data are printed.
            //  This could lead to NANs if writeEventDataToFile is called after this.
            
            ModelPtr->resetMHacceptanceParameters();
        }

        // Deprecating this - no need to write this accept data
        //if ((i % _acceptWriteFreq) == 0)
        //    writeParamAcceptRates();
    }
}


MCMC::~MCMC(void)
{
    _mcmcOutStream.close();
    _lambdaOutStream.close();
    _muOutStream.close();
    _acceptOutStream.close();
    _lambdaNodeOutStream.close();
    _eventDataOutStream.close();
}


void MCMC::setUpdateWeights(void)
{
    parWts.push_back(sttings->getUpdateRateEventNumber()); // event number
    parWts.push_back(sttings->getUpdateRateEventPosition()); // event position
    parWts.push_back(sttings->getUpdateRateEventRate()); // event rate
    parWts.push_back(sttings->getUpdateRateLambda0()); // lambda0 rate
    parWts.push_back(sttings->getUpdateRateLambdaShift()); // lambda shift
    parWts.push_back(sttings->getUpdateRateMu0()); //mu0 rate
    parWts.push_back(sttings->getUpdateRateMuShift()); // mu shift
    parWts.push_back(
        sttings->getUpdateRateNumberTimeVariablePartitions()); // freq of updates to timevarying/constant

    double sumwts = parWts[0];
    for (std::vector<double>::size_type i = 1; i < parWts.size(); i++) {
        sumwts += parWts[i];
        parWts[i] += parWts[i - 1];
    }

    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++) {
        parWts[i] /= sumwts;
    }

    //for (int i = 0; i < parWts.size(); i++)
    //  std::cout << parWts[i] << std::endl;

    // Define std::vectors to hold accept/reject data:
    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++) {
        acceptCount.push_back(0);
        rejectCount.push_back(0);
    }
}


int MCMC::pickParameterClassToUpdate(void)
{
    double rn = ranPtr->uniformRv();
    int parm = 0;
    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++) {
        if (rn <= parWts[i]) {
            parm = (int)i;
            break;
        }
    }
    return parm;
}


void MCMC::updateState(int parm)
{
    if (parm == 0) {
        ModelPtr->changeNumberOfEventsMH();
    } else if (parm == 1) {
        ModelPtr->moveEventMH();
    } else if (parm == 2) {
        ModelPtr->updateEventRateMH();
    } else if (parm == 3) {
        ModelPtr->updateLambdaInitMH();
    } else if (parm == 4) {
        ModelPtr->updateLambdaShiftMH();
    } else if (parm == 5) {
        ModelPtr->updateMuInitMH();
    } else if (parm == 6) {
        ModelPtr->updateMuShiftMH();
    } else if (parm == 7) {
        // Update time variable partition:
        std::cout << "should update isTimeVariable" << std::endl;
        ModelPtr->setAcceptLastUpdate(1);
    } else {
        // should never get here...throw exception?
        std::cout << "Bad parm to update\n" << std::endl;
    }

    if (ModelPtr->getAcceptLastUpdate() == 1) {
        acceptCount[parm]++;
    } else if ( ModelPtr->getAcceptLastUpdate() == 0 ) {
        rejectCount[parm]++;
    } else if ( ModelPtr->getAcceptLastUpdate() == -1) {
        std::cout << "failed somewhere in MH step, parm " << parm << std::endl;
        throw;
    } else {
        std::cout << "invalid accept/reject flag in model object" << std::endl;
        throw;
    }
    // reset to unmodified value
    ModelPtr->setAcceptLastUpdate(-1);
}


void MCMC::writeStateToFile()
{
    if (_firstLine) {
      writeHeaderToStream(_mcmcOutStream);
      _firstLine = false;
    } else {
      writeStateToStream(_mcmcOutStream);
    }
}


void MCMC::writeHeaderToStream(std::ostream& outStream)
{
    outStream << "generation,numevents,logprior,llbranches,eventRate,acceptRate\n";
}


void MCMC::writeStateToStream(std::ostream& outStream)
{
    outStream << ModelPtr->getGeneration()       << ","
              << ModelPtr->getNumberOfEvents()   << ","
              << ModelPtr->computeLogPrior()     << ","
              << ModelPtr->getCurrLnLBranches()  << ","
              << ModelPtr->getEventRate()        << ","
              << ModelPtr->getMHacceptanceRate() << "," << std::endl;
}


/* print state data to screen:
    current LnL
    # of events
    mean speciation rate?


*/
void MCMC::printStateData(void)
{
    std::cout << std::setw(10) << ModelPtr->getGeneration() << std::setw(10);
    std::cout << ModelPtr->getCurrLnLBranches() << std::setw(10);
    std::cout << ModelPtr->getNumberOfEvents() << std::setw(10);
    std::cout << ModelPtr->computeLogPrior() << std::setw(15);
//  std::cout << "OtherL: " << ModelPtr->computeLikelihoodBranchesByInterval() << std::setw(10);
    //std::cout << ModelPtr->getEventRate() << std::setw(10);
    std::cout << "Accp: " << ModelPtr->getMHacceptanceRate() << std::endl;
}

////////

void MCMC::writeBranchSpeciationRatesToFile(void)
{
    ModelPtr->getTreePtr()->setMeanBranchSpeciation();
    std::stringstream outdata;
    ModelPtr->getTreePtr()->writeMeanBranchSpeciationTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    _lambdaOutStream << outdata.str() << std::endl;
}


void MCMC::writeNodeSpeciationRatesToFile(void)
{
    ModelPtr->getTreePtr()->setMeanBranchSpeciation();
    std::stringstream outdata;
    ModelPtr->getTreePtr()->writeNodeSpeciationTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    _lambdaNodeOutStream << outdata.str() << std::endl;
}


void MCMC::writeBranchExtinctionRatesToFile(void)
{
    ModelPtr->getTreePtr()->setMeanBranchExtinction();
    std::stringstream outdata;
    ModelPtr->getTreePtr()->writeMeanBranchExtinctionTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    _muOutStream << outdata.str() << std::endl;
}

// TODO: Deprecate writeParamAcceptRates (not useful at present)
void MCMC::writeParamAcceptRates(void)
{
    std::ofstream& outStream = _acceptOutStream;
    for (std::vector<int>::size_type i = 0; i < acceptCount.size(); i++) {
        double rate = (double)acceptCount[i] / ((double)(acceptCount[i] +
                                                rejectCount[i]));
        outStream << rate;
        if (i < (acceptCount.size() - 1)) {
            outStream << ",";
        } else {
            outStream << "\n";
        }
    }
    for (std::vector<int>::size_type i = 0; i < acceptCount.size(); i++) {
        acceptCount[i] = 0;
        rejectCount[i] = 0;
    }
}

void MCMC::writeEventDataToFile(void)
{
    std::stringstream eventData;
    ModelPtr->getEventDataString(eventData);

    _eventDataOutStream << eventData.str() << std::endl;
}


bool MCMC::anyOutputFileExists()
{
    return fileExists(_mcmcOutFilename)       ||
           fileExists(_lambdaOutFilename)     ||
           fileExists(_muOutFilename)         ||
           fileExists(_acceptOutFilename)     ||
           fileExists(_lambdaNodeOutFilename) ||
           fileExists(_eventDataOutFilename);
}


bool MCMC::fileExists(const std::string& filename)
{
    std::ifstream inFile(filename.c_str());
    return inFile.good();
}


void MCMC::exitWithErrorOutputFileExists()
{
    std::cout << "ERROR: Analysis is set to not overwrite files.\n";
    std::cout << "Fix by removing or renaming output file(s),\n";
    std::cout << "or set \"overwrite = 1\" in the control file.\n";
    std::exit(1);
}
