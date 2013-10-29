/*
 *  mcmc.h
 *  rateshift
 *
 *  Created by Dan Rabosky on 12/2/11.
 *
 */

#ifndef MCMC_H
#define MCMC_H

#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class MbRandom;
class Model; 
class Settings;


class MCMC{

public:
	MCMC(MbRandom * ran, Model * mymodel, Settings * sp);
	~MCMC(void);
	
	void		writeStateToFile(void);
	void		printStateData(void);
	void		writeBranchSpeciationRatesToFile(void);
	void		writeBranchExtinctionRatesToFile(void);
	void		writeNodeSpeciationRatesToFile(void);
	void		writeEventDataToFile(void);
	
	int			pickParameterClassToUpdate(void);
	void		updateState(int parm);
	
	void		setUpdateWeights(void);
 
	void		writeParamAcceptRates(void);
 
	
private:
	MbRandom		*ranPtr;
	Model			*ModelPtr;
	Settings		*sttings;
	
	vector<double>	parWts;
	
	vector<int>		acceptCount;
	vector<int>		rejectCount;
	
	string			mcmcOutfile;
	string			lambdaOutfile;
	string			lambdaNodeOutfile;
	string			muOutfile;
	string			acceptFile;
	string			eventDataFile;
	
	int				_treeWriteFreq;
	int				_mcmcWriteFreq;
	int				_acceptWriteFreq;
	int				_printFreq;
	int				_NGENS;
	
};



#endif




