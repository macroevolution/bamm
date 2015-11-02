#ifndef SP_EX_DATA_WRITER_H
#define SP_EX_DATA_WRITER_H


#include "ModelDataWriter.h"
#include "SpExEventDataWriter.h"

class Settings;
class Model;


class SpExDataWriter : public ModelDataWriter
{
public:

    SpExDataWriter(Settings &settings);

    virtual void writeData(int generation, Model& model);

protected:

    SpExEventDataWriter _eventDataWriter;
};


#endif
