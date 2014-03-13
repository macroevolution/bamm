#include "TraitDataWriter.h"
#include "TraitModel.h"


TraitDataWriter::TraitDataWriter(Settings& settings, TraitModel& model) :
    DataWriter(settings, model)
{
}

void TraitDataWriter::writeSpecificEventHeaders()
{
    _eventOutputStream << ",betainit,betashift\n";
}
