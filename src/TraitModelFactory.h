#ifndef TRAIT_MODEL_FACTORY
#define TRAIT_MODEL_FACTORY


#include "ModelFactory.h"
#include "TraitModel.h"
#include "TraitDataWriter.h"

class Model;
class DataWriter;

class MbRandom;
class Settings;
class Prior;


class TraitModelFactory : public ModelFactory
{
public:

    virtual ~TraitModelFactory() {}

    virtual Model* createModel
        (MbRandom& rng, Settings& settings, Prior& prior) const;
    virtual DataWriter* createDataWriter
        (Settings& settings, Model& model) const;
};


inline Model* TraitModelFactory::createModel
    (MbRandom& rng, Settings& settings, Prior& prior) const
{
    return new TraitModel(&rng, &settings, &prior);
}


inline DataWriter* TraitModelFactory::createDataWriter
    (Settings& settings, Model& model) const
{
    return new TraitDataWriter(settings, static_cast<TraitModel&>(model));
}


#endif
