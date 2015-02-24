#include "MetropolisCoupledMCMC.h"

#include "Random.h"
#include "Settings.h"
#include "ModelFactory.h"
#include "MCMC.h"
#include "Model.h"
#include "ModelDataWriter.h"
#include "ChainSwapDataWriter.h"

#include <algorithm>
#include <thread>


MetropolisCoupledMCMC::MetropolisCoupledMCMC
    (Random& random, Settings& settings, ModelFactory* modelFactory) :
        _random(random), _settings(settings), _modelFactory(modelFactory),
        _chainSwapDataWriter(_settings)
{
    // Total number of generations to run for each chain
    _nGenerations = _settings.get<int>("numberOfGenerations");

    // MC3 settings
    _nChains = _settings.get<int>("numberOfChains");
    _deltaT = _settings.get<double>("deltaT");
    _swapPeriod = _settings.get<int>("swapPeriod");

    _coldChainIndex = 0;

    _acceptanceResetFreq = _settings.get<int>("acceptanceResetFreq");
}


MetropolisCoupledMCMC::~MetropolisCoupledMCMC()
{
    for (int i = 0; i < (int)_chains.size(); i++) {
        delete _chains[i];
    }

    delete _dataWriter;
}


void MetropolisCoupledMCMC::run()
{
     createChains();
    createDataWriter();
 
    log() << "\nRunning " << _chains.size() << " chains for "
          << _nGenerations << " generations.\n";

    log() << "\n";

    int generation = 0;
    while (generation < _nGenerations) {
        int generationEnd = std::min(generation + _swapPeriod, _nGenerations);
        runChains(generation, generationEnd);
        generation = generationEnd;
        tryChainSwap(generation);
    }
}


void MetropolisCoupledMCMC::createChains()
{
    for (int i = 0; i < _nChains; i++) {
        _chains.push_back(createMCMC(i));
    }
}


MCMC* MetropolisCoupledMCMC::createMCMC(int chainIndex) const
{
    MCMC* mcmc = new MCMC(_random, _settings, *_modelFactory);
    mcmc->model().setTemperatureMH(calculateTemperature(chainIndex, _deltaT));
    return mcmc;
}


double MetropolisCoupledMCMC::calculateTemperature(int i, double deltaT) const
{
    return 1.0 / (1.0 + deltaT * i);
}


void MetropolisCoupledMCMC::createDataWriter()
{
 
    // TODO: THere is a bug here. IF I set fossil parameters
    //  such as updateRatePreservationRate within the SpeciationExtinction
    // Parameters block, this line of code will cause an error for BAMM trait
 
    _dataWriter = _modelFactory->createModelDataWriter(_settings);
 

}


void MetropolisCoupledMCMC::runChains(int genStart, int genEnd)
{
    std::vector<std::thread> chainThreads;
    chainThreads.reserve(_chains.size());

    for (int i = 0; i < (int)_chains.size(); i++) {
        chainThreads.push_back(std::thread
            {&MetropolisCoupledMCMC::runChain, this, i, genStart, genEnd});
    }

    for (std::thread& chainThread : chainThreads) {
        chainThread.join();
    }
}


void MetropolisCoupledMCMC::runChain(int i, int genStart, int genEnd)
{
    for (int g = genStart; g < genEnd; g++) {
        _chains[i]->step();

        if (i == _coldChainIndex) {
            _dataWriter->writeData(g, _chains[i]->model());

            if (g % _acceptanceResetFreq == 0) {
                _chains[i]->model().resetMHAcceptanceParameters();
            }
        }
    }
}


void MetropolisCoupledMCMC::tryChainSwap(int generation)
{
    if ((_chains.size() == 1) || (_swapPeriod == 0) ||
            (generation % _swapPeriod != 0)) {
        return;
    }

    int chain_1, chain_2;
    chooseTwoNumbers(&chain_1, &chain_2, 0, (int)_chains.size() - 1);

    bool chainSwapAccepted = acceptChainSwap(chain_1, chain_2);

    if (chainSwapAccepted) {
        swapTemperature(chain_1, chain_2);
    }

    _chainSwapDataWriter.writeData
        (generation, _chains, chain_1, chain_2, chainSwapAccepted);
}


void MetropolisCoupledMCMC::chooseTwoNumbers(int* x, int* y, int from, int to)
{
    *x = _random.uniformInteger(from, to);

    do {
        *y = _random.uniformInteger(from, to);
    } while (*y == *x);
}


bool MetropolisCoupledMCMC::acceptChainSwap(int chain_1, int chain_2) const
{
    return _random.trueWithProbability(chainSwapProbability(chain_1, chain_2));
}


double MetropolisCoupledMCMC::chainSwapProbability
    (int chain_1, int chain_2) const
{
    Model& model_1 = _chains[chain_1]->model();
    Model& model_2 = _chains[chain_2]->model();

    double beta_1 = model_1.getTemperatureMH();
    double beta_2 = model_2.getTemperatureMH();

    double log_post_1 = calculateLogPosterior(model_1);
    double log_post_2 = calculateLogPosterior(model_2);

    double swapPosteriorRatio = std::exp
        (logSwapPosteriorRatio(beta_1, beta_2, log_post_1, log_post_2));

    return std::min(1.0, swapPosteriorRatio);
}


double MetropolisCoupledMCMC::calculateLogPosterior(Model& model) const
{
    return model.getCurrentLogLikelihood() + model.computeLogPrior();
}


double MetropolisCoupledMCMC::logSwapPosteriorRatio
    (double beta_1, double beta_2, double log_post_1, double log_post_2) const
{
    return (beta_2 - beta_1) * log_post_1 + (beta_1 - beta_2) * log_post_2;
}


void MetropolisCoupledMCMC::swapTemperature(int chain_1, int chain_2)
{
    double beta_1 = _chains[chain_1]->model().getTemperatureMH();
    double beta_2 = _chains[chain_2]->model().getTemperatureMH();

    _chains[chain_1]->model().setTemperatureMH(beta_2);
    _chains[chain_2]->model().setTemperatureMH(beta_1);

    // Properly keep track of the cold chain
    if (chain_1 == _coldChainIndex) {
        _coldChainIndex = chain_2;
    } else if (chain_2 == _coldChainIndex) {
        _coldChainIndex = chain_1;
    }
}
