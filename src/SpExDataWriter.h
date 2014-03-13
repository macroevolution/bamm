#ifndef SP_EX_DATA_WRITER_H
#define SP_EX_DATA_WRITER_H


#include "DataWriter.h"

class Settings;
class SpExModel;


class SpExDataWriter : public DataWriter
{
public:

    SpExDataWriter(Settings& settings, SpExModel& model);

private:

    void writeSpecificEventHeaders();
};


#endif
