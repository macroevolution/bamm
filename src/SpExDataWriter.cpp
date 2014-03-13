#include "SpExDataWriter.h"
#include "SpExModel.h"


SpExDataWriter::SpExDataWriter(Settings& settings, SpExModel& model) :
    DataWriter(settings, model)
{
}

void SpExDataWriter::writeSpecificEventHeaders()
{
    _eventOutputStream << ",lambdainit,lambdashift,muinit,mushift\n";
}
