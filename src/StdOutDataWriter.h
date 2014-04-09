#ifndef STD_OUT_DATA_WRITER_H
#define STD_OUT_DATA_WRITER_H


#include <string>
#include <fstream>

class Settings;
class Model;
class MCMC;


class StdOutDataWriter
{
public:

    StdOutDataWriter(Settings& settings);
    ~StdOutDataWriter();

    void writeData(int generation, MCMC& mcmc);

private:

    void writeHeader();
    std::string header();

    int _outputFreq;
    bool _headerWritten;
};


#endif
