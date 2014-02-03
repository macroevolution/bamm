/*
 *  TraitMCMC.cpp
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
  *
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "TraitMCMC.h"
#include "TraitModel.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "Settings.h"
#include "TraitBranchEvent.h"
#include "Log.h"


TraitMCMC::TraitMCMC(MbRandom* ran, TraitModel* mymodel, Settings* sp)
{
    ranPtr = ran;
    ModelPtr = mymodel;
    sttings = sp;

    //check if file exists; delete;

    // MCMC parameters - for now....//
    _mcmcOutFilename      = sttings->getMCMCoutfile();
    _betaOutFilename      = sttings->getBetaOutfile();
    _eventDataOutFilename = sttings->getEventDataOutfile();

    _writeMeanBranchLengthTrees = sttings->getWriteMeanBranchLengthTrees();

    _treeWriteFreq =    sttings->getBranchRatesWriteFreq();
    _mcmcWriteFreq =    sttings->getMCMCwriteFreq();
    _eventDataWriteFreq = sttings->getEventDataWriteFreq();
    _printFreq =        sttings->getPrintFreq();
    _NGENS =            sttings->getNGENS();

    // Open streams for writing (overwrite)
    _mcmcOutStream.open(_mcmcOutFilename.c_str());
    _eventDataOutStream.open(_eventDataOutFilename.c_str());

    if (_writeMeanBranchLengthTrees) {
        _betaOutStream.open(_betaOutFilename.c_str());
    }

    writeHeadersToOutputFiles();

    setUpdateWeights();
    ModelPtr->resetGeneration();

    for (int i = 0; i < sttings->getInitialNumberEvents(); i++)
        ModelPtr->addEventToTree();

    log() << "\nRunning MCMC chain for " << _NGENS << " generations.\n";

    log() << "\n"
          << std::setw(15) << "Generation"
          << std::setw(15) << "LogLikelihood"
          << std::setw(15) << "NumShifts"    // Why not NumEvents?
          << std::setw(15) << "LogPrior"
          << std::setw(15) << "AcceptRate"
          << "\n";

    /*run chain*/
    for (int i = 0; i < _NGENS; i++) {


        int parmToUpdate = pickParameterClassToUpdate();
        updateState(parmToUpdate); // update state



        if (_writeMeanBranchLengthTrees && (i % _treeWriteFreq == 0)) {
            writeBranchBetaRatesToFile();
        }

        if ((i % _eventDataWriteFreq) == 0){
            writeEventDataToFile();        
        }




        if ((i % _mcmcWriteFreq) == 0){
            writeStateToFile();      
        }

        


        if ((i % _printFreq == 0)){
            
            printStateData();
            ModelPtr->resetMHacceptanceParameters();
        
        }
    }

}


TraitMCMC::~TraitMCMC(void)
{
    _mcmcOutStream.close();
    _eventDataOutStream.close();

    if (_writeMeanBranchLengthTrees) {
        _betaOutStream.close();
    }
}


void TraitMCMC::setUpdateWeights(void)
{

    parWts.push_back(
        sttings->getUpdateRateEventNumber());                      // 0    event number
    parWts.push_back(
        sttings->getUpdateRateEventPosition());                    // 1    event position
    parWts.push_back(
        sttings->getUpdateRateEventRate());                        // 2    event rate
    parWts.push_back(
        sttings->getUpdateRateBeta0());                            // 3    beta0 rate
    parWts.push_back(
        sttings->getUpdateRateBetaShift());                        // 4    beta shift
    parWts.push_back(
        sttings->getUpdateRateNodeState());                        // 5    Node states


    double sumwts = parWts[0];
    for (std::vector<double>::size_type i = 1; i < parWts.size(); i++) {
        sumwts += parWts[i];
        parWts[i] += parWts[i - 1];
    }

    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++)
        parWts[i] /= sumwts;

    // Define std::vectors to hold accept/reject data:
    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++) {
        acceptCount.push_back(0);
        rejectCount.push_back(0);
    }



}


int TraitMCMC::pickParameterClassToUpdate(void)
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


void TraitMCMC::updateState(int parm)
{

    if (parm == 0)
        ModelPtr->changeNumberOfEventsMH();
    else if (parm == 1)
        ModelPtr->moveEventMH();
    else if (parm == 2)
        ModelPtr->updateEventRateMH();
    else if (parm == 3)
        ModelPtr->updateBetaMH();
    else if (parm == 4)
        ModelPtr->updateBetaShiftMH();
    else if (parm == 5)
        ModelPtr->updateNodeStateMH();
    else if (parm == 6) {
        // Update time variable partition:
        std::cout << "should update isTimeVariable" << std::endl;
        ModelPtr->setAcceptLastUpdate(1);

    } else {
        // should never get here...throw exception?
        std::cout << "Bad parm to update\n" << std::endl;
    }

    if (ModelPtr->getAcceptLastUpdate() == 1)
        acceptCount[parm]++;
    else if ( ModelPtr->getAcceptLastUpdate() == 0 )
        rejectCount[parm]++;
    else if ( ModelPtr->getAcceptLastUpdate() == -1) {
        std::cout << "failed somewhere in MH step, parm " << parm << std::endl;
        throw;
    } else {
        std::cout << "invalid accept/reject flag in model object" << std::endl;
        throw;
    }
    // reset to unmodified value
    ModelPtr->setAcceptLastUpdate(-1);

}


void TraitMCMC::writeStateToFile()
{
    writeStateToStream(_mcmcOutStream);
}


void TraitMCMC::writeHeaderToStream(std::ostream& outStream)
{
    outStream << "generation,N_shifts,logPrior,logLik,acceptRate\n";
}


void TraitMCMC::writeStateToStream(std::ostream& outStream)
{
    outStream << ModelPtr->getGeneration()     << ","
              << ModelPtr->getNumberOfEvents() << ","
              << ModelPtr->computeLogPrior()   << ","
              << ModelPtr->getCurrLnLTraits()  << ","
              << ModelPtr->getEventRate()      << ","
              << ModelPtr->getMHacceptanceRate() << std::endl;
}


/* print state data to screen:
 current LnL
 # of events
 mean speciation rate?


 */
void TraitMCMC::printStateData(void)
{
    log() << std::setw(15) << ModelPtr->getGeneration()
          << std::setw(15) << ModelPtr->getCurrLnLTraits()
          << std::setw(15) << ModelPtr->getNumberOfEvents()
          << std::setw(15) << ModelPtr->computeLogPrior()
          << std::setw(15) << ModelPtr->getMHacceptanceRate()
          << "\n";
}

////////

void TraitMCMC::writeBranchBetaRatesToFile(void)
{
    ModelPtr->getTreePtr()->setMeanBranchTraitRates();

    std::stringstream outdata;
    ModelPtr->getTreePtr()->writeMeanBranchTraitRateTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);

    outdata << ";";

    _betaOutStream << outdata.str() << std::endl;
}


void TraitMCMC::writeEventDataToFile(void)
{
    std::stringstream eventData;
    ModelPtr->getEventDataString(eventData);

    _eventDataOutStream << eventData.str() << std::endl;
}


bool TraitMCMC::anyOutputFileExists()
{
    return fileExists(_mcmcOutFilename)      ||
           fileExists(_betaOutFilename)      ||
           fileExists(_eventDataOutFilename);
}


bool TraitMCMC::fileExists(const std::string& filename)
{
    std::ifstream inFile(filename.c_str());
    return inFile.good();
}


void TraitMCMC::writeHeadersToOutputFiles()
{
    _mcmcOutStream << "generation,N_shifts,logPrior,logLik,eventRate,acceptRate\n";
    _eventDataOutStream << "generation,leftchild,rightchild,abstime," <<
        "betainit,betashift\n";
}


void TraitMCMC::exitWithErrorOutputFileExists()
{
    log(Error) << "Analysis is set to not overwrite files.\n"
               << "Fix by removing or renaming output file(s),\n"
               << "or set \"overwrite = 1\" in the control file.\n";
    std::exit(1);
}
