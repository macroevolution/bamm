#ifndef MODEL_FACTORY
#define MODEL_FACTORY


class Model;
class EventDataWriter;

class Random;
class Settings;
class Prior;


class ModelFactory
{
public:

    virtual ~ModelFactory() {}

    virtual Model* createModel(Random& random, Settings& settings) const = 0;
    virtual EventDataWriter* createEventDataWriter
        (Settings& settings) const = 0;
};


#endif
