#ifndef SP_EX_MCMC_H
#define SP_EX_MCMC_H

#include <vector>
#include <string>
#include <fstream>

#include "MCMC.h"

class MbRandom;
class Model;
class SpExModel;
class Settings;


class SpExMCMC : public MCMC
{

public:

    SpExMCMC(MbRandom* rng, Model* model, Settings* settings);
    virtual ~SpExMCMC();

private:

    virtual void setUpSpecificParameterWeights();
    virtual void updateSpecificState(int parameter);

    virtual void outputSpecificEventDataHeaders();

    virtual void outputSpecificData(int generation);
    void outputBranchSpeciationRates();
    void outputBranchExtinctionRates();

    SpExModel* _specificModel;

    std::string _lambdaOutputFileName;
    std::string _muOutputFileName;

    std::ofstream _lambdaOutputStream;
    std::ofstream _muOutputStream;

    int _treeOutputFreq;

    bool _writeMeanBranchLengthTrees;
};


#endif
