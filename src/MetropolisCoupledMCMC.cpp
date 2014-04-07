MetropolisCoupledMCMC::MetropolisCoupledMCMC(MbRandom& rng, Settings& settings,
    Prior& prior, ModelFactory* modelFactory) : _rng(rng), _settings,
    _prior(prior), _modelFactory(modelFactory)
{
    // Total number of generations to run for each chain
    _nGenerations = _settings.get<int>("numberGenerations");

    // MC3 settings
    _nChains = _settings.get<int>("nChains");
    _deltaT = _settings.get<double>("deltaT");
    _swapPeriod = _settings.get<int>("swapPeriod");

    // Track number of chain swaps proposed and accepted
    _nSwapsProposed = 0;
    _nSwapsAccepted = 0;

    _dataWriter = _modelFactory->createDataWriter(_settings);
}


MetropolisCoupledMCMC::~MetropolisCoupledMCMC()
{
    for (int i = 0; i < (int)_chains.size(); i++) {
        delete _chains[i]->model();
        delete _chains[i];
    }
}


void MetropolisCoupledMCMC::run()
{
    createChains();

    log() << "Running " << _chains.size() << " chains for "
          << _nGenerations << " generations.\n";

    for (int gen = 0; gen < _nGenerations; gen += _swapPeriod) {
        // Run each chain up to _swapPeriod steps
        for (int i = 0; i < (int)_chains.size(); i++) {
            _chains[i]->run(_swapPeriod);
        }

        if (_chains.size() > 1) {
            int chain_1, chain_2;
            chooseTwoNumbers(&chain_1, &chain_2, 0, (int)_chains.size() - 1);

            bool chainSwapAccepted = acceptChainSwap(chain_1, chain_2);

            _nSwapsProposed++;
            if (chainSwapAccepted) {
                swapTemperature(chain_1, chain_2);
                _nSwapsAccepted++;
            }
        }
    }
}


void MetropolisCoupledMCMC::createChains() const
{
    for (int i = 0; i < _nChains; i++) {
        _chains.push_back(createMCMC(i));
    }
}


MCMC* MetropolisCoupledMCMC::createMCMC(int chainIndex) const
{
    return new MCMC(_rng, createModel(chainIndex));
}


Model* MetropolisCoupledMCMC::createModel(int chainIndex) const
{
    Model* model = _modelFactory->createModel(_rng, _settings, _prior);
    model->setTemperatureMH(calculateTemperature(chainIndex, _deltaT));
    return model;
}


double MetropolisCoupledMCMC::calculateTemperature(int i, double deltaT) const
{
    return 1.0 / (1.0 + deltaT * i);
}


void MetropolisCoupledMCMC::chooseTwoNumbers(int* x, int* y, int from, int to)
{
    *x = _rng->discreteUniformRv(from, to);

    do {
        *y = _rng->discreteUniformRv(from, to);
    } while (*y == *x)
}


bool MetropolisCoupledMCMC::acceptChainSwap(int chain_1, int chain_2) const
{
    return trueWithProbability(chainSwapProbability(chain_1, chain_2));
}


bool MetropolisCoupledMCMC::trueWithProbability(double p) const
{
    return _rng->uniformRv() < p;
}


double MetropolisCoupledMCMC::chainSwapProbability
    (int chain_1, int chain_2) const
{
    Model* model_1 = _chains[chain_1]->getModel();
    Model* model_2 = _chains[chain_2]->getModel();

    double beta_1 = model_1->getTemperatureMH();
    double beta_2 = model_2->getTemperatureMH();

    double log_post_1 = calculateLogPosterior(model_1);
    double log_post_2 = calculateLogPosterior(model_2);

    double swapPosteriorRatio = std::exp
        (logSwapPosteriorRatio(beta_1, beta_2, log_post_1, log_post_2));

    return std::min(1.0, swapPosteriorRatio);
}


double MetropolisCoupledMCMC::logSwapPosteriorRatio
    (double beta_1, double beta_2, double log_post_1, double log_post_2) const
{
    return (beta_2 - beta_1) * log_post_1 + (beta_1 - beta_2) * log_post_2;
}


void MetropolisCoupledMCMC::swapTemperature(int chain_1, int chain_2)
{
    double beta_1 = _chains[chain_1]->getModel()->getTemperatureMH();
    double beta_2 = _chains[chain_2]->getModel()->getTemperatureMH();

    _chains[chain_1]->getModel()->setTemperatureMH(beta_2);
    _chains[chain_2]->getModel()->setTemperatureMH(beta_1);

    // Properly keep track of the cold chain
    if (chain_1 == _coldChainIndex) {
        _coldChainIndex = chain_2;
    } else if (chain_2 == _coldChainIndex) {
        _coldChainIndex = chain_1;
    }
}
