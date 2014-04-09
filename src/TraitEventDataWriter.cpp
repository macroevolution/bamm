#include "TraitEventDataWriter.h"
#include "EventDataWriter.h"
#include "TraitBranchEvent.h"

#include <sstream>

class Settings;
class BranchEvent;


TraitEventDataWriter::TraitEventDataWriter(Settings& settings) :
    EventDataWriter(settings)
{
}


TraitEventDataWriter::~TraitEventDataWriter()
{
}


std::string TraitEventDataWriter::eventParameters(BranchEvent* event)
{
    TraitBranchEvent* specificEvent = static_cast<TraitBranchEvent*>(event);

    std::ostringstream stringStream;
    stringStream << specificEvent->getBetaInit()   << ","
                 << specificEvent->getBetaShift();
    return stringStream.str();
}
