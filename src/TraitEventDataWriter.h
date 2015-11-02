#ifndef TRAIT_EVENT_DATA_WRITER_H
#define TRAIT_EVENT_DATA_WRITER_H


#include "EventDataWriter.h"
#include <string>

class Settings;
class BranchEvent;


class TraitEventDataWriter : public EventDataWriter
{
public:

    TraitEventDataWriter(Settings& settings);
    virtual ~TraitEventDataWriter();

private:

    virtual std::string specificHeader();
    virtual std::string eventParameters(BranchEvent* event);
};


inline std::string TraitEventDataWriter::specificHeader()
{
    return ",betainit,betashift";
}


#endif
