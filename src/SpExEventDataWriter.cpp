#include "SpExEventDataWriter.h"
#include "EventDataWriter.h"
#include "SpExBranchEvent.h"

#include <sstream>

class Settings;
class BranchEvent;


SpExEventDataWriter::SpExEventDataWriter(Settings& settings) :
    EventDataWriter(settings)
{
}


SpExEventDataWriter::~SpExEventDataWriter()
{
}


std::string SpExEventDataWriter::eventParameters(BranchEvent* event)
{
    SpExBranchEvent* specificEvent = static_cast<SpExBranchEvent*>(event);

    std::ostringstream stringStream;
    stringStream << specificEvent->getLamInit()  << ","
                 << specificEvent->getLamShift() << ","
                 << specificEvent->getMuInit()   << ","
                 << specificEvent->getMuShift();
    return stringStream.str();
}
