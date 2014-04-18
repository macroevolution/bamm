#include "MCMC.h"
#include "Random.h"
#include "Model.h"
#include "ModelFactory.h"

#include <climits>


MCMC::MCMC(Random& seeder, Settings& settings, ModelFactory& modelFactory) :
    _random(seeder.uniformInteger(0, INT_MAX))
{
    _model = modelFactory.createModel(_random, settings);
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
    if (_random.trueWithProbability(acceptanceRatio)) {
        _model->acceptProposal();
    } else {
        _model->rejectProposal();
    }
}
