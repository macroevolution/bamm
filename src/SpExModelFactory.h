#ifndef SP_EX_MODEL_FACTORY
#define SP_EX_MODEL_FACTORY


#include "ModelFactory.h"
#include "SpExModel.h"
#include "SpExDataWriter.h"

class Model;
class ModelDataWriter;

class Random;
class Settings;
class Prior;


class SpExModelFactory : public ModelFactory
{
public:

    virtual ~SpExModelFactory() {}

    virtual Model* createModel(Random& random, Settings& settings) const;
    virtual ModelDataWriter* createModelDataWriter(Settings& settings) const;
};


inline Model* SpExModelFactory::createModel
    (Random& random, Settings& settings) const
{
    return new SpExModel(random, settings);
}


inline ModelDataWriter* SpExModelFactory::createModelDataWriter
    (Settings& settings) const
{
    return new SpExDataWriter(settings);
}


#endif
