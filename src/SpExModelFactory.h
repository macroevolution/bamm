#ifndef SP_EX_MODEL_FACTORY
#define SP_EX_MODEL_FACTORY


#include "ModelFactory.h"
#include "SpExModel.h"
#include "SpExEventDataWriter.h"

class Model;
class EventDataWriter;

class MbRandom;
class Settings;
class Prior;


class SpExModelFactory : public ModelFactory
{
public:

    virtual ~SpExModelFactory() {}

    virtual Model* createModel(MbRandom& rng, Settings& settings) const;
    virtual EventDataWriter* createEventDataWriter(Settings& settings) const;
};


inline Model* SpExModelFactory::createModel
    (MbRandom& rng, Settings& settings) const
{
    return new SpExModel(&rng, &settings);
}


inline EventDataWriter* SpExModelFactory::createEventDataWriter
    (Settings& settings) const
{
    return new SpExEventDataWriter(settings);
}


#endif
