#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include <string>
#include <fstream>

class MbRandom;
class Model;
class Settings;


class MCMC
{
    typedef std::vector<double>::size_type SizeType;

public:

    MCMC(MbRandom* rng, Model* model, Settings* settings);
    virtual ~MCMC();

    void run();

protected:

    void setUpdateWeights();
    void setUpParameterWeights();
    virtual void setUpSpecificParameterWeights() = 0;

    int  chooseRandomParameter();
    void updateState(int parameter);
    virtual void updateSpecificState(int parameter) = 0;

    void outputHeaders();
    void outputMCMCHeaders();
    void outputEventDataHeaders();
    virtual void outputSpecificEventDataHeaders() = 0;
    void outputStdOutHeaders();
    void outputAcceptanceInfoHeaders();

    void outputData(int generation);
    void outputMCMCData();
    void outputEventData();
    void outputStdOutData();
    virtual void outputSpecificData(int generation) = 0;
    void outputAcceptanceInfo(int param, bool accepted);

    MbRandom* _rng;
    Model* _model;
    Settings* _settings;

    int _numGenerations;

    std::vector<double> _parameterWeights;

    std::vector<int> _acceptCount;
    std::vector<int> _rejectCount;

    std::string _mcmcOutputFileName;
    std::string _eventDataOutputFileName;
    std::string _acceptanceInfoFileName;

    std::ofstream _mcmcOutputStream;
    std::ofstream _eventDataOutputStream;
    std::ofstream _acceptanceInfoStream;

    int _mcmcOutputFreq;
    int _eventDataOutputFreq;
    int _stdOutFreq;

    bool _outputAcceptanceInfo;
};


#endif
