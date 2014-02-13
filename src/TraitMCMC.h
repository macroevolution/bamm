#ifndef TRAIT_MCMC_H
#define TRAIT_MCMC_H

#include <vector>
#include <string>
#include <fstream>

#include "MCMC.h"

class MbRandom;
class Model;
class TraitModel;
class Settings;


class TraitMCMC : public MCMC
{

public:

    TraitMCMC(MbRandom* rng, Model* model, Settings* settings);
    virtual ~TraitMCMC();

private:

    virtual void setUpSpecificParameterWeights();
    virtual void updateSpecificState(int parameter);

    virtual void outputSpecificEventDataHeaders();

    virtual void outputSpecificData(int generation);
    virtual void outputBranchBetaRates();

    TraitModel* _specificModel;

    bool _writeMeanBranchLengthTrees;
    int _treeOutputFreq;

    std::string _betaOutputFileName;
    std::ofstream _betaOutputStream;
};


#endif
