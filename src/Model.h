#ifndef MODEL_H
#define MODEL_H


class MbRandom;
class Tree;
class Settings;
class Prior;


class Model
{

public:

    static double mhColdness;

    Model(MbRandom* rng, Tree* tree, Settings* settings, Prior* prior);
    virtual ~Model();

    Tree* getTreePtr();

    void incrementGeneration();
    int getGeneration();
    void resetGeneration();

    void setMoveSizeScale(double x);
    void setUpdateEventRateScale(double x);

    double getPoissonRatePrior();
    void setPoissonRatePrior(double x);

    double getEventRate();
    void setEventRate(double x);

    int getAcceptLastUpdate();
    void setAcceptLastUpdate(int x);

protected:

    MbRandom* _rng;
    Tree* _tree;
    Settings* _settings;
    Prior* _prior;

    int _gen;

    // Parameters for MCMC proposals
    double _scale;    // scale for moving event
    double _updateEventRateScale;
    double _localGlobalMoveRatio;

    double _poissonRatePrior;
    double _eventRate;    // Poisson rate

    int _acceptCount;
    int _rejectCount;
    int _acceptLast;    // true if last generation was accept; false otherwise
    // 0 = last was rejected; 1 = accepted; -1 = not set.

    double _lastDeletedEventMapTime;    // map time of last deleted event
};


inline Tree* Model::getTreePtr()
{
    return _tree;
}


inline void Model::incrementGeneration()
{
    _gen++;
}


inline int Model::getGeneration()
{
    return _gen;
}


inline void Model::resetGeneration()
{
    _gen = 0;    // to be used after TraitPreBurnin
}


inline void Model::setMoveSizeScale(double x)
{
    _scale = x;
}


inline void Model::setUpdateEventRateScale(double x)
{
    _updateEventRateScale = x;
}


inline double Model::getPoissonRatePrior()
{
    return _poissonRatePrior;
}


inline void Model::setPoissonRatePrior(double x)
{
    _poissonRatePrior = x;
}


inline double Model::getEventRate()
{
    return _eventRate;
}


inline void Model::setEventRate(double x)
{
    _eventRate = x;
}


inline int Model::getAcceptLastUpdate()
{
    return _acceptLast;
}


inline void Model::setAcceptLastUpdate(int x)
{
    _acceptLast = x;
}


#endif
