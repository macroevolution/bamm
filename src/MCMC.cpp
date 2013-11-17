#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

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
    mcmcOutfile             =   sttings->getMCMCoutfile();
    lambdaOutfile           =   sttings->getLambdaOutfile();
    muOutfile               =   sttings->getMuOutfile();
    acceptFile              =   sttings->getAcceptrateOutfile();
    lambdaNodeOutfile       =   sttings->getLambdaNodeOutfile();
    eventDataFile           =   sttings->getEventDataOutfile();

    _treeWriteFreq =        sttings->getTreeWriteFreq();
    _eventDataWriteFreq =   sttings->getEventDataWriteFreq();
    _mcmcWriteFreq =        sttings->getMCMCwriteFreq();
    _acceptWriteFreq =      sttings->getAcceptWriteFreq();
    _printFreq =            sttings->getPrintFreq();
    _NGENS =                sttings->getNGENS();

    _firstLine = true; // Print header for the first line of output


    std::ifstream outStream(mcmcOutfile.c_str());
    if (outStream) {
        std::cout << "Output file for MCMC data exists: " << std::endl;
        std::cout << std::setw(30) << " overwriting " << mcmcOutfile << std::endl;
        outStream.close();

        std::string filedelete("rm ");
        filedelete.append(mcmcOutfile);

        system(filedelete.c_str());

    }

    //check if file exists; delete;
    std::ifstream outStream2(lambdaOutfile.c_str());
    if (outStream2) {
        std::cout << "Output file for lambda exists: " << std::endl;
        std::cout << std::setw(30) << " overwriting " << lambdaOutfile << std::endl;
        outStream2.close();

        std::string filedelete("rm ");
        filedelete.append(lambdaOutfile);

        system(filedelete.c_str());

    }

    std::ifstream outStream3(muOutfile.c_str());
    if (outStream3) {
        std::cout << "Output file for mu exists: " << std::endl;
        std::cout << std::setw(30) << " overwriting " << muOutfile << std::endl;
        outStream3.close();

        std::string filedelete("rm ");
        filedelete.append(muOutfile);

        system(filedelete.c_str());

    }

    std::ifstream outStream4(acceptFile.c_str());
    if (outStream4) {
        std::cout << "Output file for acceptrates exists: " << std::endl;
        std::cout << std::setw(30) << " overwriting " << acceptFile << std::endl;
        outStream4.close();
        std::string filedelete("rm ");
        filedelete.append(acceptFile);
        system(filedelete.c_str());
    }

    std::ifstream outStream5(lambdaNodeOutfile.c_str());
    if (outStream5) {
        std::cout << "Output file for lambdaNodeData: " << std::endl;
        std::cout << std::setw(30) << " overwriting " << lambdaNodeOutfile << std::endl;
        outStream5.close();
        std::string filedelete("rm ");
        filedelete.append(lambdaNodeOutfile);
        system(filedelete.c_str());
    }

    std::ifstream outStream6(eventDataFile.c_str());
    if (outStream6) {
        std::cout << "Output file for event data: " << std::endl;
        std::cout << std::setw(30) << " overwriting " << eventDataFile << std::endl;
        outStream6.close();
        std::string filedelete("rm ");
        filedelete.append(eventDataFile);
        system(filedelete.c_str());
    }


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

            writeBranchSpeciationRatesToFile();
            writeBranchExtinctionRatesToFile();
            writeNodeSpeciationRatesToFile();
            //writeEventDataToFile();
            //mymodel->getTreePtr()->printNodeBranchRates();
        }

        if ((i % _eventDataWriteFreq) == 0)

            writeEventDataToFile();


        if ((i % _mcmcWriteFreq) == 0)
            writeStateToFile();



        if ((i % _printFreq == 0))
            printStateData();

        if ((i % _acceptWriteFreq) == 0)
            writeParamAcceptRates();

    }

}


MCMC::~MCMC(void)
{

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

    if (parm == 0)
        ModelPtr->changeNumberOfEventsMH();
    else if (parm == 1)
        ModelPtr->moveEventMH();
    else if (parm == 2)
        ModelPtr->updateEventRateMH();
    else if (parm == 3)
        ModelPtr->updateLambdaInitMH();
    else if (parm == 4)
        ModelPtr->updateLambdaShiftMH();
    else if (parm == 5)
        ModelPtr->updateMuInitMH();
    else if (parm == 6)
        ModelPtr->updateMuShiftMH();
    else if (parm == 7) {
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


void MCMC::writeStateToFile()
{
    // TODO: Move opening to constructor, closing to destructor
    std::ofstream outStream(mcmcOutfile.c_str(), std::ofstream::app);

    if (_firstLine) {
      writeHeaderToStream(outStream);
      _firstLine = false;
    }
    else
      writeStateToStream(outStream);

    outStream.close();
}


void MCMC::writeHeaderToStream(std::ostream& outStream)
{
    outStream << "generation,numevents,logprior,llbranches,eventRate\n";
}


void MCMC::writeStateToStream(std::ostream& outStream)
{
    outStream << ModelPtr->getGeneration()      << ","
              << ModelPtr->getNumberOfEvents()  << ","
              << ModelPtr->computeLogPrior()    << ","
              << ModelPtr->getCurrLnLBranches() << ","
			  << ModelPtr->getEventRate() << std::endl;
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

    std::string outname = lambdaOutfile;

    ModelPtr->getTreePtr()->setMeanBranchSpeciation();

    std::stringstream outdata;

    ModelPtr->getTreePtr()->writeMeanBranchSpeciationTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);

    outdata << ";";

    std::ofstream outStream;
    outStream.open(outname.c_str(), std::ofstream::app);
    outStream << outdata.str() << std::endl;
    outStream.close();

}


void MCMC::writeNodeSpeciationRatesToFile(void)
{

    std::string outname = lambdaNodeOutfile;

    ModelPtr->getTreePtr()->setMeanBranchSpeciation();

    std::stringstream outdata;

    ModelPtr->getTreePtr()->writeNodeSpeciationTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);

    outdata << ";";

    std::ofstream outStream;
    outStream.open(outname.c_str(), std::ofstream::app);
    outStream << outdata.str() << std::endl;
    outStream.close();

}


void MCMC::writeBranchExtinctionRatesToFile(void)
{

    std::string outname = muOutfile;


    ModelPtr->getTreePtr()->setMeanBranchExtinction();
    std::stringstream outdata;
    ModelPtr->getTreePtr()->writeMeanBranchExtinctionTree(
        ModelPtr->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    std::ofstream outStream;
    outStream.open(outname.c_str(), std::ofstream::app);
    outStream << outdata.str() << std::endl;
    outStream.close();

}

void MCMC::writeParamAcceptRates(void)
{

    std::ofstream outStream;
    outStream.open(acceptFile.c_str(), std::ofstream::app);
    for (std::vector<int>::size_type i = 0; i < acceptCount.size(); i++) {
        double rate = (double)acceptCount[i] / ((double)(acceptCount[i] +
                                                rejectCount[i]));
        outStream << rate;
        if (i < (acceptCount.size() - 1))
            outStream << ",";
        else
            outStream << "\n";

    }
    outStream.close();

    for (std::vector<int>::size_type i = 0; i < acceptCount.size(); i++) {
        acceptCount[i] = 0;
        rejectCount[i] = 0;
    }

}

void MCMC::writeEventDataToFile(void)
{

    std::string outname = eventDataFile;

    std::stringstream eventData;
    ModelPtr->getEventDataString(eventData);

    std::ofstream outstream;
    outstream.open(eventDataFile.c_str(), std::ofstream::app);
    outstream << eventData.str() << std::endl;
    outstream.close();

}
