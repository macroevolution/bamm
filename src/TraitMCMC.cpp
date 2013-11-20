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


TraitMCMC::TraitMCMC(MbRandom* ran, TraitModel* mymodel, Settings* sp)
{

    std::cout << "Initializing Trait MCMC object..." << std::endl;

    ranPtr = ran;
    ModelPtr = mymodel;
    sttings = sp;

    //check if file exists; delete;

    // MCMC parameters - for now....//
    _mcmcOutFilename      = sttings->getMCMCoutfile();
    _betaOutFilename      = sttings->getBetaOutfile();
    _nodeStateOutFilename = sttings->getNodeStateOutfile();
    _acceptOutFilename    = sttings->getAcceptrateOutfile();
    _eventDataOutFilename = sttings->getEventDataOutfile();

    _treeWriteFreq =    sttings->getBranchRatesWriteFreq();
    _mcmcWriteFreq =    sttings->getMCMCwriteFreq();
    _eventDataWriteFreq = sttings->getEventDataWriteFreq();
    _acceptWriteFreq =  sttings->getAcceptWriteFreq();
    _printFreq =        sttings->getPrintFreq();
    _NGENS =            sttings->getNGENS();

    bool fileOverwrite = sttings->getOverwrite();

    if (!fileOverwrite) {
        if (anyOutputFileExists()) {
            exitWithErrorOutputFileExists();
        }
    }

    // Open streams for writing (overwrite)
    _mcmcOutStream.open(_mcmcOutFilename.c_str());
    _betaOutStream.open(_betaOutFilename.c_str());
    _nodeStateOutStream.open(_nodeStateOutFilename.c_str());
    _acceptOutStream.open(_acceptOutFilename.c_str());
    _eventDataOutStream.open(_eventDataOutFilename.c_str());

    writeHeadersToOutputFiles();

    setUpdateWeights();
    ModelPtr->resetGeneration();

    for (int i = 0; i < sttings->getInitialNumberEvents(); i++)
        ModelPtr->addEventToTree();

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

            writeBranchBetaRatesToFile();
            writeNodeStatesToFile();

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

        

        
        //if ((i % _acceptWriteFreq) == 0)
        //  writeParamAcceptRates();

    }

}


TraitMCMC::~TraitMCMC(void)
{

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

    parWts.push_back(
        sttings->getUpdateRateNumberTimeVariablePartitions());     // 6 freq of updates to timevarying/constant

    double sumwts = parWts[0];
    for (std::vector<double>::size_type i = 1; i < parWts.size(); i++) {
        sumwts += parWts[i];
        parWts[i] += parWts[i - 1];
    }

    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++)
        parWts[i] /= sumwts;

    //for (int i = 0; i < parWts.size(); i++)
    //  std::cout << parWts[i] << std::endl;

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
    outStream << "generation,numevents,logprior,lltraits,acceptRate\n";
}


void TraitMCMC::writeStateToStream(std::ostream& outStream)
{
    outStream << ModelPtr->getGeneration()     << ","
              << ModelPtr->getNumberOfEvents() << ","
              << ModelPtr->computeLogPrior()   << ","
              << ModelPtr->getCurrLnLTraits()  << ","
              << ModelPtr->getMHacceptanceRate() << std::endl;
}


/* print state data to screen:
 current LnL
 # of events
 mean speciation rate?


 */
void TraitMCMC::printStateData(void)
{

    std::cout << std::setw(10) << ModelPtr->getGeneration() << std::setw(10);
    std::cout << ModelPtr->getCurrLnLTraits() << std::setw(10);
    std::cout << ModelPtr->getNumberOfEvents() << std::setw(10);
    std::cout << ModelPtr->computeLogPrior() << std::setw(15);
    //  std::cout << "OtherL: " << ModelPtr->computeLikelihoodBranchesByInterval() << std::setw(10);
    //std::cout << ModelPtr->getEventRate() << std::setw(10);
    std::cout << "Accp: " << ModelPtr->getMHacceptanceRate() << std::endl;

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


void TraitMCMC::writeNodeStatesToFile(void)
{
    std::stringstream outdata;

    ModelPtr->getTreePtr()->writeBranchPhenotypes
        (ModelPtr->getTreePtr()->getRoot(), outdata);

    outdata << ";";

    _nodeStateOutStream << outdata.str() << std::endl;
}

void TraitMCMC::writeParamAcceptRates(void)
{
    std::ofstream& outStream = _acceptOutStream;
    for (std::vector<int>::size_type i = 0; i < acceptCount.size(); i++) {
        double rate = (double)acceptCount[i] / ((double)(acceptCount[i] +
                                                rejectCount[i]));
        outStream << rate;
        if (i < (acceptCount.size() - 1))
            outStream << ",";
        else
            outStream << "\n";

    }

    for (std::vector<int>::size_type i = 0; i < acceptCount.size(); i++) {
        acceptCount[i] = 0;
        rejectCount[i] = 0;
    }
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
           fileExists(_nodeStateOutFilename) ||
           fileExists(_acceptOutFilename)    ||
           fileExists(_eventDataOutFilename);
}


bool TraitMCMC::fileExists(const std::string& filename)
{
    std::ifstream inFile(filename.c_str());
    return inFile.good();
}


void TraitMCMC::writeHeadersToOutputFiles()
{
    _mcmcOutStream << "generation,numevents,logprior,lltraits,acceptRate\n";
    _eventDataOutStream << "generation,leftchild,rightchild,abstime," <<
        "betainit,betashift\n";
}


void TraitMCMC::exitWithErrorOutputFileExists()
{
    std::cout << "ERROR: Analysis is set to not overwrite files.\n";
    std::cout << "Fix by removing or renaming output file(s),\n";
    std::cout << "or set \"overwrite = 1\" in the control file.\n";
    std::exit(1);
}
