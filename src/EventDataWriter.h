#ifndef EVENT_DATA_WRITER_H
#define EVENT_DATA_WRITER_H

#include <iostream>
#include <fstream>

class Settings;
class Model;
class BranchEvent;
class Node;


class EventDataWriter
{
public:

    EventDataWriter(Settings& settings);
    virtual ~EventDataWriter();

    void writeData(int generation, Model& model);

protected:

    void writeHeaderOnce();
    void writeHeader();
    std::string header();
    virtual std::string specificHeader() = 0;

    void writeEventData(int generation, Model& model);
    void writeRootEvent(int generation, Model& model);
    void writeEvents(int generation, Model& model);

    void writeEvent(int generation, BranchEvent* event);
    std::string leftNodeName(BranchEvent* event);
    std::string rightNodeName(BranchEvent* event);
    virtual std::string eventParameters(BranchEvent* event) = 0;

    std::string _outputFileName;
    std::ofstream _outputStream;
    int _outputFreq;

    bool _headerWritten;
};


#endif
