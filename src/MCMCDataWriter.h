#ifndef MCMC_DATA_WRITER_H
#define MCMC_DATA_WRITER_H


#include <string>
#include <fstream>

class Settings;
class Model;


class MCMCDataWriter
{
public:

    MCMCDataWriter(Settings& settings);
    ~MCMCDataWriter();

    void writeData(int generation, Model& model);

private:

    void initializeStream();
    void writeHeader();
    std::string header();

    std::string _outputFileName;
    int _outputFreq;
    
    std::ofstream _outputStream;

    bool _hasPreservationRate;

};


#endif
