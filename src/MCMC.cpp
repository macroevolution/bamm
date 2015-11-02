#include "MCMC.h"
#include "Random.h"
#include "Model.h"
#include "ModelFactory.h"

#include <climits>


// Choose a random number up to INT_MAX - 1, not INT_MAX,
// because MbRandom adds 1 internally, causing an overflow
MCMC::MCMC(Random& seeder, Settings& settings, ModelFactory& modelFactory) :
    _random(seeder.uniformInteger(0, INT_MAX - 1))
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
    //std::cout << _model->getCurrentLogLikelihood() << "\tActual: " << _model->computeLogLikelihood() << std::endl;
    
    //double logL = _model->computeLogLikelihood();
  
     _model->proposeNewState();

    double acceptanceRatio = _model->acceptanceRatio();
    if (_random.trueWithProbability(acceptanceRatio)) {
        _model->acceptProposal();
    } else {
        _model->rejectProposal();
    }
}
