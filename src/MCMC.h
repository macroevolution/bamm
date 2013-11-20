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

    bool anyOutputFileExists();
    bool fileExists(const std::string& filename);
    void writeHeadersToOutputFiles();
    void exitWithErrorOutputFileExists();

    MbRandom* ranPtr;
    Model*    ModelPtr;
    Settings* sttings;

    std::vector<double> parWts;

    std::vector<int> acceptCount;
    std::vector<int> rejectCount;

    std::string _mcmcOutFilename;
    std::string _lambdaOutFilename;
    std::string _lambdaNodeOutFilename;
    std::string _muOutFilename;
    std::string _acceptOutFilename;
    std::string _eventDataOutFilename;

    std::ofstream _mcmcOutStream;
    std::ofstream _lambdaOutStream;
    std::ofstream _muOutStream;
    std::ofstream _acceptOutStream;
    std::ofstream _lambdaNodeOutStream;
    std::ofstream _eventDataOutStream;

    int _treeWriteFreq;
    int _eventDataWriteFreq;
    int _mcmcWriteFreq;
    int _acceptWriteFreq;
    int _printFreq;
    int _NGENS;
};


#endif
