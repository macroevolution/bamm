/*
 *  TraitMCMC.h
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.

 *
 */

#ifndef TRAIT_MCMC_H
#define TRAIT_MCMC_H

#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

class MbRandom;
class TraitModel;
class Settings;


class TraitMCMC
{

public:

    TraitMCMC(MbRandom* ran, TraitModel* mymodel, Settings* sp);
    ~TraitMCMC(void);

    void writeStateToFile(void);
    void printStateData(void);
    void writeBranchBetaRatesToFile(void);
    void writeNodeStatesToFile(void);
    void writeEventDataToFile(void);

    int  pickParameterClassToUpdate(void);
    void updateState(int parm);

    void setUpdateWeights(void);

    void writeParamAcceptRates(void);

private:

    MbRandom* ranPtr;
    TraitModel* ModelPtr;
    Settings* sttings;

    vector<double> parWts;

    vector<int> acceptCount;
    vector<int> rejectCount;

    string mcmcOutfile;
    string betaOutfile;
    string nodeStateOutfile;
    string acceptFile;
    string eventDataFile;

    int _treeWriteFreq;
    int _mcmcWriteFreq;
    int _eventDataWriteFreq;
    int _acceptWriteFreq;
    int _printFreq;
    int _NGENS;
};


#endif
