



#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "MCMC.h"
#include "model.h"
#include "MbRandom.h"
#include "node.h"
#include "Settings.h"
 
MCMC::MCMC(MbRandom * ran, Model * mymodel, Settings * sp){
	
	cout << "Initializing MCMC object..." << endl;
	
	ranPtr = ran;
	ModelPtr = mymodel;
	sttings = sp;
	
	//check if file exists; delete;
	
	// MCMC parameters - for now....// 
	mcmcOutfile				=	sttings->getMCMCoutfile();
	lambdaOutfile			=	sttings->getLambdaOutfile();
	muOutfile				=	sttings->getMuOutfile();
	acceptFile				=	sttings->getAcceptrateOutfile();
	lambdaNodeOutfile		=	sttings->getLambdaNodeOutfile();
	eventDataFile			=	sttings->getEventDataOutfile();
	
	_treeWriteFreq =		sttings->getTreeWriteFreq();
	_eventDataWriteFreq = 	sttings->getEventDataWriteFreq();
	_mcmcWriteFreq =		sttings->getMCMCwriteFreq();
	_acceptWriteFreq =		sttings->getAcceptWriteFreq();
	_printFreq =			sttings->getPrintFreq();
	_NGENS =				sttings->getNGENS();
	
	
	ifstream outStream(mcmcOutfile.c_str());
	if (outStream){
		cout << "Output file for MCMC data exists: " << endl;
		cout << setw(30) << " overwriting " << mcmcOutfile << endl;
 		outStream.close();
		
		string filedelete("rm ");
		filedelete.append(mcmcOutfile);
		
		system(filedelete.c_str());
		
	}
	
	//check if file exists; delete;
	ifstream outStream2(lambdaOutfile.c_str());
	if (outStream2){
		cout << "Output file for lambda exists: " << endl;
		cout << setw(30) << " overwriting " << lambdaOutfile << endl;	
 		outStream2.close();
		
		string filedelete("rm ");
		filedelete.append(lambdaOutfile);
		
		system(filedelete.c_str());
		
	}	
	
	ifstream outStream3(muOutfile.c_str());
	if (outStream3){
		cout << "Output file for mu exists: " << endl;
		cout << setw(30) << " overwriting " << muOutfile << endl;
 		outStream3.close();
		
		string filedelete("rm ");
		filedelete.append(muOutfile);
		
		system(filedelete.c_str());
		
	}		
	
	ifstream outStream4(acceptFile.c_str());
	if (outStream4){
		cout << "Output file for acceptrates exists: " << endl;
		cout << setw(30) << " overwriting " << acceptFile << endl;
		outStream4.close();
		string filedelete("rm ");
		filedelete.append(acceptFile);
		system(filedelete.c_str());
	}
	
	ifstream outStream5(lambdaNodeOutfile.c_str());
	if (outStream5){
		cout << "Output file for lambdaNodeData: " << endl;
		cout << setw(30) << " overwriting " << lambdaNodeOutfile << endl;
		outStream5.close();
		string filedelete("rm ");
		filedelete.append(lambdaNodeOutfile);
		system(filedelete.c_str());
	}	
	
	ifstream outStream6(eventDataFile.c_str());
	if (outStream6){
		cout << "Output file for event data: " << endl;
		cout << setw(30) << " overwriting " << eventDataFile << endl;
		outStream6.close();
		string filedelete("rm ");
		filedelete.append(eventDataFile);
		system(filedelete.c_str());	
	}
	
	
	setUpdateWeights();
	ModelPtr->resetGeneration();
	
	for (int i = 0; i < sttings->getInitialNumberEvents(); i++)
		ModelPtr->addEventToTree();
	
	cout << "MCMC object successfully initialized." << endl << endl;
	cout << "Running MCMC chain for " << _NGENS << " generations." << endl << endl;
	
	cout << setw(10) << "Generation" << setw(10) << "lnLik"<<  setw(10);
	cout << "N_shifts" << setw(15) << "LogPrior" << setw(15) << "acceptRate" << endl;
 	
	/*run chain*/
	for (int i = 0; i < _NGENS; i++){

 
		int parmToUpdate = pickParameterClassToUpdate();
		updateState(parmToUpdate); // update state
		
		

		if ((i % _treeWriteFreq) == 0){
			
			writeBranchSpeciationRatesToFile();
			writeBranchExtinctionRatesToFile();
			writeNodeSpeciationRatesToFile();
			//writeEventDataToFile();		
			//mymodel->getTreePtr()->printNodeBranchRates();
		}
		
		if ((i % _eventDataWriteFreq) == 0){
			
			writeEventDataToFile();
			
		}
		
		if ((i % _mcmcWriteFreq) == 0){
			writeStateToFile();

		}
		
		
		if ((i % _printFreq == 0)){
			printStateData();
		}
		
		if ((i % _acceptWriteFreq) == 0){
			writeParamAcceptRates();
		}
 
	}

}


MCMC::~MCMC(void){

}


void MCMC::setUpdateWeights(void){

	parWts.push_back(sttings->getUpdateRateEventNumber()); // event number
	parWts.push_back(sttings->getUpdateRateEventPosition()); // event position
	parWts.push_back(sttings->getUpdateRateEventRate()); // event rate 
	parWts.push_back(sttings->getUpdateRateLambda0()); // lambda0 rate
	parWts.push_back(sttings->getUpdateRateLambdaShift()); // lambda shift
	parWts.push_back(sttings->getUpdateRateMu0()); //mu0 rate
	parWts.push_back(sttings->getUpdateRateMuShift()); // mu shift
	parWts.push_back(sttings->getUpdateRateNumberTimeVariablePartitions()); // freq of updates to timevarying/constant

	double sumwts = parWts[0];
	for (vector<double>::size_type i = 1; i < parWts.size(); i++){
		sumwts += parWts[i];
		parWts[i] += parWts[i-1];	
	}

	for (vector<double>::size_type i = 0; i < parWts.size(); i++)
		parWts[i] /= sumwts;
	
	//for (int i = 0; i < parWts.size(); i++)
	//	cout << parWts[i] << endl;
	
	// Define vectors to hold accept/reject data:
	for (vector<double>::size_type i = 0; i < parWts.size(); i++){
		acceptCount.push_back(0);
		rejectCount.push_back(0);
	}
	
	
	
}


int MCMC::pickParameterClassToUpdate(void){
	double rn = ranPtr->uniformRv();
	int parm = 0;
	for (vector<double>::size_type i = 0; i < parWts.size(); i++){
		if (rn <= parWts[i]){
			parm = i;
			break;
		}
	}
	return parm;
}


void MCMC::updateState(int parm){
	
	if (parm == 0){
		ModelPtr->changeNumberOfEventsMH();
	}else if (parm == 1){
		ModelPtr->moveEventMH();
	}else if (parm == 2){
		ModelPtr->updateEventRateMH();
 	}else if (parm == 3){
		ModelPtr->updateLambdaInitMH();
	}else if(parm == 4){
		ModelPtr->updateLambdaShiftMH();
	}else if (parm == 5){
		ModelPtr->updateMuInitMH();
	}else if(parm == 6){
		ModelPtr->updateMuShiftMH();
	}else if(parm == 7){
		// Update time variable partition:
		cout << "should update isTimeVariable" << endl;
		ModelPtr->setAcceptLastUpdate(1);
		
	} else{
		// should never get here...throw exception?
		cout << "Bad parm to update\n" << endl;
	}
		
	if (ModelPtr->getAcceptLastUpdate() == 1){
		acceptCount[parm]++;
	}else if ( ModelPtr->getAcceptLastUpdate() == 0 ){
		rejectCount[parm]++;
	}else if ( ModelPtr->getAcceptLastUpdate() == -1){
		cout << "failed somewhere in MH step, parm " << parm << endl;
		throw;
	}else{
		cout << "invalid accept/reject flag in model object" << endl;
		throw;
	}
	// reset to unmodified value
	ModelPtr->setAcceptLastUpdate(-1);
	
}



void MCMC::writeStateToFile(void){
	
	
	//replace file...
 
	
	ofstream outStream;
	outStream.open(mcmcOutfile.c_str(), ofstream::app);
	outStream << ModelPtr->getGeneration() << ","  << ModelPtr->getNumberOfEvents() << ",";
	outStream << ModelPtr->computeLogPrior() << ",";
	outStream << ModelPtr->getCurrLnLBranches() << ",";
	outStream << ModelPtr->getEventRate() << "," << ModelPtr->getRootEvent()->getLamInit() << ",";
	outStream << ModelPtr->getRootEvent()->getLamShift() << ",";
	outStream << ModelPtr->getRootEvent()->getMuInit() << ",";
	outStream << ModelPtr->getRootEvent()->getMuShift() << endl;
	outStream.close();
 
}


/* print state data to screen:
	current LnL
	# of events
	mean speciation rate?
 
 
*/
void MCMC::printStateData(void){
	
	cout << setw(10) << ModelPtr->getGeneration() << setw(10);
	cout << ModelPtr->getCurrLnLBranches() << setw(10);
	cout << ModelPtr->getNumberOfEvents() << setw(10);
	cout << ModelPtr->computeLogPrior() << setw(15);
//	cout << "OtherL: " << ModelPtr->computeLikelihoodBranchesByInterval() << setw(10);
	//cout << ModelPtr->getEventRate() << setw(10);
	cout << "Accp: " << ModelPtr->getMHacceptanceRate() << endl;
 
}

////////

void MCMC::writeBranchSpeciationRatesToFile(void){
	
	string outname = lambdaOutfile;
	
	ModelPtr->getTreePtr()->setMeanBranchSpeciation();
	
	stringstream outdata;
 
	ModelPtr->getTreePtr()->writeMeanBranchSpeciationTree(ModelPtr->getTreePtr()->getRoot(), outdata);

	outdata << ";";
	
	ofstream outStream;
	outStream.open(outname.c_str(), ofstream::app);
	outStream << outdata.str() << endl;
	outStream.close();
	
}


void MCMC::writeNodeSpeciationRatesToFile(void){
	
	string outname = lambdaNodeOutfile;
	
	ModelPtr->getTreePtr()->setMeanBranchSpeciation();
	
	stringstream outdata;
	
	ModelPtr->getTreePtr()->writeNodeSpeciationTree(ModelPtr->getTreePtr()->getRoot(), outdata);
	
	outdata << ";";
	
	ofstream outStream;
	outStream.open(outname.c_str(), ofstream::app);
	outStream << outdata.str() << endl;
	outStream.close();
	
}


void MCMC::writeBranchExtinctionRatesToFile(void){
	
	string outname = muOutfile;
	
	
	ModelPtr->getTreePtr()->setMeanBranchExtinction();
	stringstream outdata;
	ModelPtr->getTreePtr()->writeMeanBranchExtinctionTree(ModelPtr->getTreePtr()->getRoot(), outdata);
	outdata << ";";
	
	ofstream outStream;
	outStream.open(outname.c_str(), ofstream::app);
	outStream << outdata.str() << endl;
	outStream.close();
	
}

void MCMC::writeParamAcceptRates(void){

	ofstream outStream;
	outStream.open(acceptFile.c_str(), ofstream::app);
	for (vector<int>::size_type i = 0; i < acceptCount.size(); i++){
		double rate = (double)acceptCount[i] / ((double)(acceptCount[i] + rejectCount[i]));
		outStream << rate;
		if (i < (acceptCount.size() - 1)){
			outStream << ",";		
		}else{
			outStream << "\n";
		}

	}
	outStream.close();
	
	for (vector<int>::size_type i = 0; i < acceptCount.size(); i++){
		acceptCount[i] = 0;
		rejectCount[i] = 0;
	}	

}

void MCMC::writeEventDataToFile(void){
	
	string outname = eventDataFile;
	
	stringstream eventData;
	ModelPtr->getEventDataString(eventData);
	
	ofstream outstream;
	outstream.open(eventDataFile.c_str(), ofstream::app);
	outstream << eventData.str() << endl;
	outstream.close();
	
}












