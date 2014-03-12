#include "SpExMCMC.h"
#include "MCMC.h"
#include "MbRandom.h"
#include "Model.h"
#include "SpExModel.h"
#include "Settings.h"
#include "Tree.h"
#include "Log.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>


SpExMCMC::SpExMCMC(MbRandom* rng, Model* model, Settings* settings) :
    MCMC(rng, model, settings)
{
    _specificModel = static_cast<SpExModel*>(model);

    _lambdaOutputFileName = _settings->getLambdaOutfile();
    _muOutputFileName     = _settings->getMuOutfile();

    // TODO: Rename to _branchRatesOutputFreq
    _treeOutputFreq = _settings->getBranchRatesWriteFreq();

    // Determine whether to write mean branch length trees
    _writeMeanBranchLengthTrees = _settings->getWriteMeanBranchLengthTrees();
    if (_writeMeanBranchLengthTrees) {
        _lambdaOutputStream.open(_lambdaOutputFileName.c_str());
        _muOutputStream.open(_muOutputFileName.c_str());
    }
}


SpExMCMC::~SpExMCMC()
{
    if (_writeMeanBranchLengthTrees) {
        _lambdaOutputStream.close();
        _muOutputStream.close();
    }
}


void SpExMCMC::outputSpecificEventDataHeaders()
{
    _eventDataOutputStream << ",lambdainit,lambdashift,muinit,mushift\n";
}


void SpExMCMC::outputSpecificData(int generation)
{
    if (_writeMeanBranchLengthTrees && (generation % _treeOutputFreq == 0)) {
        outputBranchSpeciationRates();
        outputBranchExtinctionRates();
    }
}


void SpExMCMC::outputBranchSpeciationRates()
{
    _specificModel->getTreePtr()->setMeanBranchSpeciation();

    std::stringstream outdata;
    _specificModel->getTreePtr()->writeMeanBranchSpeciationTree
        (_specificModel->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    _lambdaOutputStream << outdata.str() << std::endl;
}


void SpExMCMC::outputBranchExtinctionRates()
{
    _specificModel->getTreePtr()->setMeanBranchExtinction();

    std::stringstream outdata;
    _specificModel->getTreePtr()->writeMeanBranchExtinctionTree
        (_specificModel->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    _muOutputStream << outdata.str() << std::endl;
}
