#ifndef MODEL_FACTORY
#define MODEL_FACTORY


class Model;
class DataWriter;

class MbRandom;
class Settings;
class Prior;


class ModelFactory
{
public:

    virtual ~ModelFactory() {}

    virtual Model* createModel
        (MbRandom& rng, Settings& settings, Prior& prior) const = 0;
    virtual DataWriter* createDataWriter
        (Settings& settings, Model& model) const = 0;
};


#endif
