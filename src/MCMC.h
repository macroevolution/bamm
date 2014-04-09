#ifndef MCMC_H
#define MCMC_H


class MbRandom;
class Settings;
class Model;
class ModelFactory;


class MCMC
{
public:

    MCMC(MbRandom& rng, Settings& settings, ModelFactory& modelFactory);
    ~MCMC();

    void run(int generations);
    void step();

    Model& model();

protected:

    MbRandom& _rng;
    Model* _model;
};


inline Model& MCMC::model()
{
    return *_model;
}


#endif
