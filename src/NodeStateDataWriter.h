#ifndef NODE_STATE_DATA_WRITER_H
#define NODE_STATE_DATA_WRITER_H


#include <string>
#include <fstream>

class Settings;
class TraitModel;


class NodeStateDataWriter
{
public:

    NodeStateDataWriter(Settings& settings);
    ~NodeStateDataWriter();

    void writeData(int generation, TraitModel& model);

private:

    void initializeStream();

    std::string _outputFileName;
    int _outputFreq;

    std::ofstream _outputStream;
};


#endif
