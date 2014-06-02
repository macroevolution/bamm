#ifndef SP_EX_MODEL_FACTORY
#define SP_EX_MODEL_FACTORY


#include "ModelFactory.h"
#include "SpExModel.h"
#include "SpExEventDataWriter.h"

class Model;
class EventDataWriter;

class Random;
class Settings;
class Prior;


class SpExModelFactory : public ModelFactory
{
public:

    virtual ~SpExModelFactory() {}

    virtual Model* createModel(Random& random, Settings& settings) const;
    virtual EventDataWriter* createEventDataWriter(Settings& settings) const;
};


inline Model* SpExModelFactory::createModel
    (Random& random, Settings& settings) const
{
    return new SpExModel(random, settings);
}


inline EventDataWriter* SpExModelFactory::createEventDataWriter
    (Settings& settings) const
{
    return new SpExEventDataWriter(settings);
}


#endif
