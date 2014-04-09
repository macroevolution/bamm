#ifndef MCMC_DATA_WRITER_H
#define MCMC_DATA_WRITER_H


#include <string>
#include <fstream>

class Settings;
class Model;
class MCMC;


class MCMCDataWriter
{
public:

    MCMCDataWriter(Settings& settings);
    ~MCMCDataWriter();

    void writeData(int generation, MCMC& mcmc);

private:

    void initializeStream();
    void writeHeader();
    std::string header();

    std::string _outputFileName;
    int _outputFreq;

    std::ofstream _outputStream;
};


#endif
