#include "TraitDataWriter.h"
#include "ModelDataWriter.h"
#include "NodeStateDataWriter.h"
#include "TraitModel.h"

class Settings;
class Model;


TraitDataWriter::TraitDataWriter(Settings &settings) :
    ModelDataWriter(settings), _eventDataWriter(settings),
    _nodeStateDataWriter(settings)
{
}


void TraitDataWriter::writeData(int generation, Model& model)
{
    ModelDataWriter::writeData(generation, model);
    _eventDataWriter.writeData(generation, model);
    _nodeStateDataWriter.writeData(generation, static_cast<TraitModel&>(model));
}
