#ifndef ACCEPTANCE_DATA_WRITER_H
#define ACCEPTANCE_DATA_WRITER_H


#include <string>
#include <fstream>

class Settings;
class Model;


class AcceptanceDataWriter
{
public:

    AcceptanceDataWriter(const Settings& settings);
    ~AcceptanceDataWriter();

    void writeData(Model& model);

private:

    void initializeStream();
    void writeHeader();
    std::string header();

    bool _shouldOutputData;

    std::string _outputFileName;
    std::ofstream _outputStream;
};


#endif
