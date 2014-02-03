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
#include <iosfwd>

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
    void writeEventDataToFile(void);

    int  pickParameterClassToUpdate(void);
    void updateState(int parm);

    void setUpdateWeights(void);

private:

    void writeHeaderToStream(std::ostream& outStream);
    void writeStateToStream(std::ostream& outStream);

    bool anyOutputFileExists();
    bool fileExists(const std::string& filename);
    void writeHeadersToOutputFiles();
    void exitWithErrorOutputFileExists();

    MbRandom* ranPtr;
    TraitModel* ModelPtr;
    Settings* sttings;

    std::vector<double> parWts;

    std::vector<int> acceptCount;
    std::vector<int> rejectCount;

    std::string _mcmcOutFilename;
    std::string _betaOutFilename;
    std::string _eventDataOutFilename;

    std::ofstream _mcmcOutStream;
    std::ofstream _betaOutStream;
    std::ofstream _eventDataOutStream;

    bool _writeMeanBranchLengthTrees;

    int _treeWriteFreq;
    int _mcmcWriteFreq;
    int _eventDataWriteFreq;
    int _printFreq;
    int _NGENS;
};


#endif
