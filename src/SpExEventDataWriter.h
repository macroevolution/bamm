#ifndef SP_EX_EVENT_DATA_WRITER_H
#define SP_EX_EVENT_DATA_WRITER_H


#include "EventDataWriter.h"
#include <string>

class Settings;
class BranchEvent;


class SpExEventDataWriter : public EventDataWriter
{
public:

    SpExEventDataWriter(Settings& settings);
    virtual ~SpExEventDataWriter();

private:

    virtual std::string specificHeader();
    virtual std::string eventParameters(BranchEvent* event);
};


inline std::string SpExEventDataWriter::specificHeader()
{
    return ",lambdainit,lambdashift,muinit,mushift";
}


#endif
