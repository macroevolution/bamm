#ifndef MODEL_FACTORY
#define MODEL_FACTORY


class Model;
class EventDataWriter;
class ModelDataWriter;

class Random;
class Settings;
class Prior;


class ModelFactory
{
public:

    virtual ~ModelFactory() {}

    virtual Model* createModel(Random& random, Settings& settings) const = 0;
    virtual ModelDataWriter* createModelDataWriter
        (Settings& settings) const = 0;
};


#endif
