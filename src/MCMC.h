/*
 *  MCMC.h
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
#include <iosfwd>

class MbRandom;
class Model;
class Settings;


class MCMC
{

public:

    MCMC(MbRandom* ran, Model* mymodel, Settings* sp);
    ~MCMC();

    void writeStateToFile();
    void printStateData();
    void writeBranchSpeciationRatesToFile();
    void writeBranchExtinctionRatesToFile();
    void writeNodeSpeciationRatesToFile();
    void writeEventDataToFile();

    int  pickParameterClassToUpdate();
    void updateState(int parm);

    void setUpdateWeights();

    void writeParamAcceptRates();

private:

    void writeHeaderToStream(std::ostream& outStream);
    void writeStateToStream(std::ostream& outStream);

    MbRandom* ranPtr;
    Model*    ModelPtr;
    Settings* sttings;

    std::vector<double> parWts;

    std::vector<int> acceptCount;
    std::vector<int> rejectCount;

    std::string mcmcOutfile;
    std::string lambdaOutfile;
    std::string lambdaNodeOutfile;
    std::string muOutfile;
    std::string acceptFile;
    std::string eventDataFile;

    int _treeWriteFreq;
    int _eventDataWriteFreq;
    int _mcmcWriteFreq;
    int _acceptWriteFreq;
    int _printFreq;
    int _NGENS;

    bool _firstLine;
};


#endif
