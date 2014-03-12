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
#include <cstdlib>


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
