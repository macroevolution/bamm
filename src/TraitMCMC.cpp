#include "TraitMCMC.h"
#include "MCMC.h"
#include "MbRandom.h"
#include "Model.h"
#include "TraitModel.h"
#include "Settings.h"
#include "Node.h"
#include "Tree.h"
#include "Log.h"

#include <vector>
#include <iostream>
#include <sstream>


TraitMCMC::TraitMCMC(MbRandom* rng, Model* model, Settings* settings) :
    MCMC(rng, model, settings)
{
    _specificModel = static_cast<TraitModel*>(model);

    _betaOutputFileName = _settings->getBetaOutfile();

    _treeOutputFreq = _settings->getBranchRatesWriteFreq();

    _writeMeanBranchLengthTrees = _settings->getWriteMeanBranchLengthTrees();
    if (_writeMeanBranchLengthTrees) {
        _betaOutputStream.open(_betaOutputFileName.c_str());
    }
}


TraitMCMC::~TraitMCMC(void)
{
    if (_writeMeanBranchLengthTrees) {
        _betaOutputStream.close();
    }
}


void TraitMCMC::setUpSpecificParameterWeights()
{
    _parameterWeights.push_back(_settings->getUpdateRateBeta0());
    _parameterWeights.push_back(_settings->getUpdateRateBetaShift());
    _parameterWeights.push_back(_settings->getUpdateRateNodeState());
    // TODO: Was node state depricated?
}


void TraitMCMC::updateSpecificState(int parameter)
{
    if (parameter == 3)
        _specificModel->updateBetaMH();
    else if (parameter == 4)
        _specificModel->updateBetaShiftMH();
    else if (parameter == 5)
        _specificModel->updateNodeStateMH();
    else if (parameter == 6) {
        // Update time variable partition
        log() << "Should update isTimeVariable\n";
        _specificModel->setAcceptLastUpdate(1);
    } else {
        // Should never get here
        log(Error) << "Bad parameter to update\n";
        std::exit(1);
    }
}


void TraitMCMC::outputSpecificEventDataHeaders()
{
    _eventDataOutputStream << ",betainit,betashift\n";
}


void TraitMCMC::outputSpecificData(int generation)
{
    if (_writeMeanBranchLengthTrees && (generation % _treeOutputFreq == 0)) {
        outputBranchBetaRates();
    }
}


void TraitMCMC::outputBranchBetaRates()
{
    _specificModel->getTreePtr()->setMeanBranchTraitRates();

    std::stringstream outdata;
    _specificModel->getTreePtr()->writeMeanBranchTraitRateTree
        (_specificModel->getTreePtr()->getRoot(), outdata);
    outdata << ";";

    _betaOutputStream << outdata.str() << std::endl;
}
