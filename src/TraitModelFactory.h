#ifndef TRAIT_MODEL_FACTORY
#define TRAIT_MODEL_FACTORY


#include "ModelFactory.h"
#include "TraitModel.h"
#include "TraitEventDataWriter.h"

class Model;
class EventDataWriter;

class Random;
class Settings;
class Prior;


class TraitModelFactory : public ModelFactory
{
public:

    virtual ~TraitModelFactory() {}

    virtual Model* createModel(Random& random, Settings& settings) const;
    virtual EventDataWriter* createEventDataWriter(Settings& settings) const;
};


inline Model* TraitModelFactory::createModel
    (Random& random, Settings& settings) const
{
    return new TraitModel(random, &settings);
}


inline EventDataWriter* TraitModelFactory::createEventDataWriter
    (Settings& settings) const
{
    return new TraitEventDataWriter(settings);
}


#endif
