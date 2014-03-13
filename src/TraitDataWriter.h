#ifndef TRAIT_DATA_WRITER_H
#define TRAIT_DATA_WRITER_H


#include "DataWriter.h"

class Settings;
class TraitModel;


class TraitDataWriter : public DataWriter
{
public:

    TraitDataWriter(Settings& settings, TraitModel& model);

private:

    void writeSpecificEventHeaders();
};


#endif
