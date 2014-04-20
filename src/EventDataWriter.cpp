#include "EventDataWriter.h"
#include "Settings.h"
#include "Model.h"
#include "BranchEvent.h"
#include "Node.h"

#include <iostream>


EventDataWriter::EventDataWriter(Settings& settings) :
    _outputFileName(settings.get("eventDataOutfile")),
    _outputFreq(settings.get<int>("eventDataWriteFreq")),
    _headerWritten(false)
{
    if (_outputFreq > 0) {
        _outputStream.open(_outputFileName.c_str());
    }
}



EventDataWriter::~EventDataWriter()
{
    if (_outputFreq > 0) {
        _outputStream.close();
    }
}


void EventDataWriter::writeData(int generation, Model& model)
{
    if (_outputFreq == 0 || generation % _outputFreq != 0) {
        return;
    }

    writeHeaderOnce();
    writeEventData(generation, model);
}


void EventDataWriter::writeHeaderOnce()
{
    if (_headerWritten) {
        return;
    }

    writeHeader();
    _headerWritten = true;
}


void EventDataWriter::writeHeader()
{
    _outputStream << header() << specificHeader() << std::endl;
}


std::string EventDataWriter::header()
{
    return "generation,leftchild,rightchild,abstime";
}


void EventDataWriter::writeEventData(int generation, Model& model)
{
    writeRootEvent(generation, model);
    writeEvents(generation, model);
}


void EventDataWriter::writeRootEvent(int generation, Model& model)
{
    writeEvent(generation, model.getRootEvent());
}


void EventDataWriter::writeEvents(int generation, Model& model)
{
    EventSet& events = model.events();

    EventSet::iterator it;
    for (it = events.begin(); it != events.end(); ++it) {
        writeEvent(generation, *it);
    }
}


void EventDataWriter::writeEvent(int generation, BranchEvent* event)
{
    _outputStream << generation               << ","
                  << leftNodeName(event)      << ","
                  << rightNodeName(event)     << ","
                  << event->getAbsoluteTime() << ","
                  << eventParameters(event)   << std::endl;
}


std::string EventDataWriter::leftNodeName(BranchEvent* event)
{
    Node* eventNode = event->getEventNode();
    if (eventNode->getIsTip()) {
        return eventNode->getName();
    } else {
        return eventNode->getRandomLeftTipNode()->getName();
    }
}


std::string EventDataWriter::rightNodeName(BranchEvent* event)
{
    Node* eventNode = event->getEventNode();
    if (eventNode->getIsTip()) {
        return "NA";
    } else {
        return eventNode->getRandomRightTipNode()->getName();
    }
}
