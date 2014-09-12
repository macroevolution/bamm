#include "ModelDataWriter.h"

#include "Settings.h"
#include "StdOutDataWriter.h"
#include "MCMCDataWriter.h"
#include "AcceptanceDataWriter.h"

class Model;


ModelDataWriter::ModelDataWriter(Settings &settings) :
    _settings(settings), _stdOutDataWriter(_settings),
    _mcmcDataWriter(_settings), _acceptanceDataWriter(_settings)
{
}


ModelDataWriter::~ModelDataWriter()
{
}


void ModelDataWriter::writeData(int generation, Model& model)
{
    _stdOutDataWriter.writeData(generation, model);
    _mcmcDataWriter.writeData(generation, model);
    _acceptanceDataWriter.writeData(model);
}
