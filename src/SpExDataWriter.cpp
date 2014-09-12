#include "SpExDataWriter.h"
#include "ModelDataWriter.h"

class Settings;
class Model;


SpExDataWriter::SpExDataWriter(Settings &settings) :
    ModelDataWriter(settings), _eventDataWriter(settings)
{
}


void SpExDataWriter::writeData(int generation, Model& model)
{
    ModelDataWriter::writeData(generation, model);
    _eventDataWriter.writeData(generation, model);
}
