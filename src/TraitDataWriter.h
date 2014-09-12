#ifndef TRAIT_DATA_WRITER_H
#define TRAIT_DATA_WRITER_H


#include "ModelDataWriter.h"
#include "TraitEventDataWriter.h"
#include "NodeStateDataWriter.h"

class Settings;
class Model;


class TraitDataWriter : public ModelDataWriter
{
public:

    TraitDataWriter(Settings &settings);

    virtual void writeData(int generation, Model& model);

protected:

    TraitEventDataWriter _eventDataWriter;
    NodeStateDataWriter _nodeStateDataWriter;
};


#endif
