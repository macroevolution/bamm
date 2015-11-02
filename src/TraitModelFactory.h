#ifndef TRAIT_MODEL_FACTORY
#define TRAIT_MODEL_FACTORY


#include "ModelFactory.h"
#include "TraitModel.h"
#include "TraitDataWriter.h"

class Model;
class ModelDataWriter;

class Random;
class Settings;
class Prior;


class TraitModelFactory : public ModelFactory
{
public:

    virtual ~TraitModelFactory() {}

    virtual Model* createModel(Random& random, Settings& settings) const;
    virtual ModelDataWriter* createModelDataWriter(Settings& settings) const;
};


inline Model* TraitModelFactory::createModel
    (Random& random, Settings& settings) const
{
    return new TraitModel(random, settings);
}


inline ModelDataWriter* TraitModelFactory::createModelDataWriter
    (Settings& settings) const
{
     return new TraitDataWriter(settings);
}


#endif
