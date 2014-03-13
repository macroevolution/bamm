#ifndef MCMC_H
#define MCMC_H


#include <vector>
#include <string>
#include <fstream>

class Model;
class DataWriter;


class MCMC
{
public:

    MCMC(Model& model, int numberOfGenerations, DataWriter& dataWriter);
    ~MCMC();

    void run();

protected:

    Model& _model;
    int _numberOfGenerations;
    DataWriter& _dataWriter;
};


#endif
