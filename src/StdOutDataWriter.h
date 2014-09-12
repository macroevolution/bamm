#ifndef STD_OUT_DATA_WRITER_H
#define STD_OUT_DATA_WRITER_H


#include <string>
#include <fstream>

class Settings;
class Model;


class StdOutDataWriter
{
public:

    StdOutDataWriter(Settings& settings);
    ~StdOutDataWriter();

    void writeData(int generation, Model& model);

private:

    void writeHeader();
    std::string header();

    int _outputFreq;
    bool _headerWritten;
};


#endif
