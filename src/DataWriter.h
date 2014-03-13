#ifndef DATA_WRITER_H
#define DATA_WRITER_H


#include <string>
#include <fstream>

class Settings;
class Model;


class DataWriter
{
public:

    DataWriter(Settings& settings, Model& model);
    virtual ~DataWriter();

    void writeData(int generation);

protected:

    void initializeStreams();

    void writeHeaders();
    void writeMCMCHeaders();
    void writeEventHeaders();
    virtual void writeSpecificEventHeaders() = 0;
    void writeStdoutHeaders();
    void writeAcceptanceHeaders();

    void writeMCMCData(int generation);
    void writeEventData(int generation);
    void writeStdoutData(int generation);
    void writeAcceptanceData();

    Model& _model;

    bool _streamsInitialized;
    bool _headersWritten;

    std::string _mcmcOutputFileName;
    std::string _eventOutputFileName;
    std::string _acceptanceOutputFileName;
    bool _outputAcceptanceInfo;

    int _mcmcOutputFreq;
    int _eventOutputFreq;
    int _stdoutOutputFreq;

    std::ofstream _mcmcOutputStream;
    std::ofstream _eventOutputStream;
    std::ofstream _acceptanceOutputStream;
};


#endif
