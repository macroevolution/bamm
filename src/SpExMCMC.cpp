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


void SpExMCMC::setUpSpecificParameterWeights()
{
    _parameterWeights.push_back(_settings->getUpdateRateLambda0());
    _parameterWeights.push_back(_settings->getUpdateRateLambdaShift());
    _parameterWeights.push_back(_settings->getUpdateRateMu0());
    _parameterWeights.push_back(_settings->getUpdateRateMuShift());
}


void SpExMCMC::updateSpecificState(int parameter)
{
    if (parameter == 3) {
        _specificModel->updateLambdaInitMH();
    } else if (parameter == 4) {
        _specificModel->updateLambdaShiftMH();
    } else if (parameter == 5) {
        _specificModel->updateMuInitMH();
    } else if (parameter == 6) {
        _specificModel->updateMuShiftMH();
    } else if (parameter == 7) {
        // Update time variable partition (TODO: Is this being used?)
        log() << "Should update isTimeVariable\n";
        _specificModel->setAcceptLastUpdate(1);
    } else {
        // Should never get here
        log(Error) << "Bad parameter to update\n";
        std::exit(1);
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
