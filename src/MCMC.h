#ifndef MCMC_H
#define MCMC_H


class MbRandom;
class Model;


class MCMC
{
public:

    MCMC(MbRandom& rng, Model& model);
    ~MCMC();

    void run(int generations);
    Model& model();

protected:

    MbRandom& _rng;
    Model& _model;
};


inline Model& MCMC::model()
{
    return _model;
}


#endif
