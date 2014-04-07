#include "MCMC.h"
#include "MbRandom.h"
#include "Model.h"


MCMC::MCMC(MbRandom& rng, Model& model) : _rng(rng), _model(model)
{
}


MCMC::~MCMC()
{
}


void MCMC::run(int generations)
{
    for (int g = 0; g < generations; g++) {
        _model.proposeNewState();

        double acceptanceRatio = _model.acceptanceRatio();
        if (_rng.uniformRv() < acceptanceRatio) {
            _model.acceptProposal();
        } else {
            _model.rejectProposal();
        }
    }
}
