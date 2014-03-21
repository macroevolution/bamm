#ifndef MCMC_H
#define MCMC_H


#include <vector>
#include <string>
#include <fstream>

class MbRandom;
class Model;
class DataWriter;


class MCMC
{
public:

    MCMC(MbRandom& rng, Model& model, int numberOfGenerations,
        DataWriter& dataWriter);
    ~MCMC();

    void run();

protected:

    MbRandom& _rng;
    Model& _model;

    int _numberOfGenerations;
    DataWriter& _dataWriter;
};


#endif
