#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "MCMC.h"
#include "SpExModel.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "Settings.h"
#include "BranchEvent.h"
#include "Log.h"


MCMC::MCMC(MbRandom* ran, SpExModel* mymodel, Settings* sp)
{
    ranPtr = ran;
    ModelPtr = mymodel;
    sttings = sp;

    //check if file exists; delete;

    // MCMC parameters - for now....//
    _mcmcOutFilename       = sttings->getMCMCoutfile();
    _lambdaOutFilename     = sttings->getLambdaOutfile();
    _muOutFilename         = sttings->getMuOutfile();
    _eventDataOutFilename  = sttings->getEventDataOutfile();

    _writeMeanBranchLengthTrees = sttings->getWriteMeanBranchLengthTrees();

    _treeWriteFreq =        sttings->getBranchRatesWriteFreq();
    _eventDataWriteFreq =   sttings->getEventDataWriteFreq();
    _mcmcWriteFreq =        sttings->getMCMCwriteFreq();
    _printFreq =            sttings->getPrintFreq();
    _NGENS =                sttings->getNGENS();

    // Open streams for writing
    _mcmcOutStream.open(_mcmcOutFilename.c_str());
    _eventDataOutStream.open(_eventDataOutFilename.c_str());

    if (_writeMeanBranchLengthTrees) {
        _lambdaOutStream.open(_lambdaOutFilename.c_str());
        _muOutStream.open(_muOutFilename.c_str());
    }

    writeHeadersToOutputFiles();

    setUpdateWeights();
    ModelPtr->resetGeneration();

    for (int i = 0; i < sttings->getInitialNumberEvents(); i++) {
        ModelPtr->addEventToTree();
    }

    log() << "\nRunning MCMC chain for " << _NGENS << " generations.\n";

    log() << "\n"
          << std::setw(15) << "Generation"
          << std::setw(15) << "LogLikelihood"
          << std::setw(15) << "NumShifts"    // Why not NumEvents?
          << std::setw(15) << "LogPrior"
          << std::setw(15) << "AcceptRate"
          << std::endl;

    /*run chain*/
    for (int i = 0; i < _NGENS; i++) {
        int parmToUpdate = pickParameterClassToUpdate();
        updateState(parmToUpdate); // update state

        if (_writeMeanBranchLengthTrees && (i % _treeWriteFreq == 0)) {
            writeBranchSpeciationRatesToFile();
            writeBranchExtinctionRatesToFile();
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
    }
}


MCMC::~MCMC(void)
{
    _mcmcOutStream.close();
    _eventDataOutStream.close();

    if (_writeMeanBranchLengthTrees) {
        _lambdaOutStream.close();
        _muOutStream.close();
    }
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
    writeStateToStream(_mcmcOutStream);
}


void MCMC::writeStateToStream(std::ostream& outStream)
{
    outStream << ModelPtr->getGeneration()       << ","
              << ModelPtr->getNumberOfEvents()   << ","
              << ModelPtr->computeLogPrior()     << ","
              << ModelPtr->getCurrentLogLikelihood()  << ","
              << ModelPtr->getEventRate()        << ","
              << ModelPtr->getMHacceptanceRate() << std::endl;
}


/* print state data to screen:
    current LnL
    # of events
    mean speciation rate?


*/
void MCMC::printStateData(void)
{
    log() << std::setw(15) << ModelPtr->getGeneration()
          << std::setw(15) << ModelPtr->getCurrentLogLikelihood()
          << std::setw(15) << ModelPtr->getNumberOfEvents()
          << std::setw(15) << ModelPtr->computeLogPrior()
          << std::setw(15) << ModelPtr->getMHacceptanceRate()
          << std::endl;
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


void MCMC::writeBranchExtinctionRatesToFile(void)
{
    ModelPtr->getTreePtr()->setMeanBranchExtinction();
    std::stringstream outdata;
    ModelPtr->getTreePtr()->writeMeanBranchExtinctionTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    _muOutStream << outdata.str() << std::endl;
}


void MCMC::writeEventDataToFile(void)
{
    std::stringstream eventData;
    ModelPtr->getEventDataString(eventData);

    _eventDataOutStream << eventData.str() << std::endl;
}


void MCMC::writeHeadersToOutputFiles()
{
    _mcmcOutStream << "generation,N_shifts,logPrior,logLik," <<
        "eventRate,acceptRate" << std::endl;
    _eventDataOutStream << "generation,leftchild,rightchild,abstime," <<
        "lambdainit,lambdashift,muinit,mushift" << std::endl;
}
