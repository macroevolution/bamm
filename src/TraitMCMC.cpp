/*
 *  TraitMCMC.cpp
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
  *
 */

#include "TraitMCMC.h"

 

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "TraitMCMC.h"
#include "TraitModel.h"
#include "MbRandom.h"
#include "node.h"
#include "Settings.h"

TraitMCMC::TraitMCMC(MbRandom * ran, TraitModel * mymodel, Settings * sp){
	
	cout << "Initializing Trait MCMC object..." << endl;
	
	ranPtr = ran;
	ModelPtr = mymodel;
	sttings = sp;
	
	//check if file exists; delete;
	
	// MCMC parameters - for now....// 
	mcmcOutfile				=	sttings->getMCMCoutfile();
	betaOutfile				=	sttings->getBetaOutfile();
	nodeStateOutfile		=	sttings->getNodeStateOutfile();
	acceptFile				=	sttings->getAcceptrateOutfile();
 	eventDataFile			=	sttings->getEventDataOutfile();
	
	_treeWriteFreq =	sttings->getTreeWriteFreq();
	_mcmcWriteFreq =	sttings->getMCMCwriteFreq();
	_eventDataWriteFreq = sttings->getEventDataWriteFreq();
	_acceptWriteFreq =	sttings->getAcceptWriteFreq();
	_printFreq =		sttings->getPrintFreq();
	_NGENS =			sttings->getNGENS();
	
	
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
	ifstream outStream2(betaOutfile.c_str());
	if (outStream2){
		cout << "Output file for beta exists: " << endl;
		cout << setw(30) << " overwriting " << betaOutfile << endl;	
 		outStream2.close();
		
		string filedelete("rm ");
		filedelete.append(betaOutfile);
		
		system(filedelete.c_str());
		
	}	
	
	ifstream outStream3(nodeStateOutfile.c_str());
	if (outStream3){
		cout << "Output file for node states exists: " << endl;
		cout << setw(30) << " overwriting " << nodeStateOutfile << endl;
 		outStream3.close();
		
		string filedelete("rm ");
		filedelete.append(nodeStateOutfile);
		
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
		}
		
		if ((i % _acceptWriteFreq) == 0){
			writeParamAcceptRates();
		}
		
	}
	
}


TraitMCMC::~TraitMCMC(void){
	
}


void TraitMCMC::setUpdateWeights(void){
	
	parWts.push_back(sttings->getUpdateRateEventNumber());						// 0	event number
	parWts.push_back(sttings->getUpdateRateEventPosition());					// 1	event position
	parWts.push_back(sttings->getUpdateRateEventRate());						// 2	event rate 
	parWts.push_back(sttings->getUpdateRateBeta0());							// 3	beta0 rate
	parWts.push_back(sttings->getUpdateRateBetaShift());						// 4	beta shift
	parWts.push_back(sttings->getUpdateRateNodeState());						// 5	Node states
 
	parWts.push_back(sttings->getUpdateRateNumberTimeVariablePartitions());		// 6 freq of updates to timevarying/constant
	
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


int TraitMCMC::pickParameterClassToUpdate(void){
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


void TraitMCMC::updateState(int parm){
	
	if (parm == 0){
		ModelPtr->changeNumberOfEventsMH();
	}else if (parm == 1){
		ModelPtr->moveEventMH();
	}else if (parm == 2){
		ModelPtr->updateEventRateMH();
	}else if (parm == 3){
		ModelPtr->updateBetaMH();
	}else if(parm == 4){
		ModelPtr->updateBetaShiftMH();
	}else if (parm == 5){
		ModelPtr->updateNodeStateMH();
	}else if(parm == 6){
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



void TraitMCMC::writeStateToFile(void){
	
	
	//replace file...
	
	
	ofstream outStream;
	outStream.open(mcmcOutfile.c_str(), ofstream::app);
	outStream << ModelPtr->getGeneration() << ","  << ModelPtr->getNumberOfEvents() << ",";
	outStream << ModelPtr->computeLogPrior() << ",";
	outStream << ModelPtr->getCurrLnLTraits() << ",";
	outStream << ModelPtr->getEventRate() << "," << ModelPtr->getRootEvent()->getBetaInit() << ",";
	outStream << ModelPtr->getRootEvent()->getBetaShift() << ",";
	outStream << ModelPtr->getRootEvent()->getEventNode()->getTraitValue() << endl;
	outStream.close();
	
}


/* print state data to screen:
 current LnL
 # of events
 mean speciation rate?
 
 
 */
void TraitMCMC::printStateData(void){
	
	cout << setw(10) << ModelPtr->getGeneration() << setw(10);
	cout << ModelPtr->getCurrLnLTraits() << setw(10);
	cout << ModelPtr->getNumberOfEvents() << setw(10);
	cout << ModelPtr->computeLogPrior() << setw(15);
	//	cout << "OtherL: " << ModelPtr->computeLikelihoodBranchesByInterval() << setw(10);
	//cout << ModelPtr->getEventRate() << setw(10);
	cout << "Accp: " << ModelPtr->getMHacceptanceRate() << endl;
	
}

////////

void TraitMCMC::writeBranchBetaRatesToFile(void){
	
	string outname = betaOutfile;
	
	ModelPtr->getTreePtr()->setMeanBranchTraitRates();
	
	stringstream outdata;
	ModelPtr->getTreePtr()->writeMeanBranchTraitRateTree(ModelPtr->getTreePtr()->getRoot(), outdata);
 
	outdata << ";";
	
	ofstream outStream;
	outStream.open(outname.c_str(), ofstream::app);
	outStream << outdata.str() << endl;
	outStream.close();
	
}


void TraitMCMC::writeNodeStatesToFile(void){
	
	string outname =  nodeStateOutfile;

	stringstream outdata;
	
	ModelPtr->getTreePtr()->writeBranchPhenotypes(ModelPtr->getTreePtr()->getRoot(), outdata);
	
	outdata << ";";
	
	ofstream outStream;
	outStream.open(outname.c_str(), ofstream::app);
	outStream << outdata.str() << endl;
	outStream.close();
	
}

void TraitMCMC::writeParamAcceptRates(void){
	
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

void TraitMCMC::writeEventDataToFile(void){
	
	string outname = eventDataFile;
	
	stringstream eventData;
	ModelPtr->getEventDataString(eventData);
	
	ofstream outstream;
	outstream.open(eventDataFile.c_str(), ofstream::app);
	outstream << eventData.str() << endl;
	outstream.close();
	
}












