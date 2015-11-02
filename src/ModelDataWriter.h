#ifndef MODEL_DATA_WRITER_H
#define MODEL_DATA_WRITER_H


#include "StdOutDataWriter.h"
#include "MCMCDataWriter.h"
#include "AcceptanceDataWriter.h"

class Settings;
class Model;


class ModelDataWriter
{
public:

    ModelDataWriter(Settings &settings);
    virtual ~ModelDataWriter();

    virtual void writeData(int generation, Model& model);

protected:

    Settings &_settings;

    StdOutDataWriter _stdOutDataWriter;
    MCMCDataWriter _mcmcDataWriter;
    AcceptanceDataWriter _acceptanceDataWriter;
};


#endif
