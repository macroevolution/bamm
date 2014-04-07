#ifndef SP_EX_MODEL_FACTORY
#define SP_EX_MODEL_FACTORY


#include "ModelFactory.h"
#include "SpExModel.h"
#include "SpExDataWriter.h"

class Model;
class DataWriter;

class MbRandom;
class Settings;
class Prior;


class SpExModelFactory : public ModelFactory
{
public:

    virtual ~SpExModelFactory() {}

    virtual Model* createModel
        (MbRandom& rng, Settings& settings, Prior& prior) const;
    virtual DataWriter* createDataWriter
        (Settings& settings, Model& model) const;
};


inline Model* SpExModelFactory::createModel
    (MbRandom& rng, Settings& settings, Prior& prior) const
{
    return new SpExModel(&rng, &settings, &prior);
}


inline DataWriter* SpExModelFactory::createDataWriter
    (Settings& settings, Model& model) const
{
    return new SpExDataWriter(settings, static_cast<SpExModel&>(model));
}


#endif
