#include "MCMC.h"
#include "MbRandom.h"
#include "Model.h"
#include "ModelFactory.h"


MCMC::MCMC(MbRandom& rng, Settings& settings, ModelFactory& modelFactory) :
    _rng(rng)
{
    _model = modelFactory.createModel(_rng, settings);
}


MCMC::~MCMC()
{
    delete _model;
}


void MCMC::run(int generations)
{
    for (int g = 0; g < generations; g++) {
        step();
    }
}


void MCMC::step()
{
    _model->proposeNewState();

    double acceptanceRatio = _model->acceptanceRatio();
    if (_rng.uniformRv() < acceptanceRatio) {
        _model->acceptProposal();
    } else {
        _model->rejectProposal();
    }
}
