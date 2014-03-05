#ifndef MODEL_H
#define MODEL_H


#include "BranchEvent.h"

#include <set>
#include <iosfwd>

class MbRandom;
class Tree;
class Settings;
class Prior;
class Node;


class Model
{

public:

    typedef std::set<BranchEvent*, BranchEvent::PtrCompare> EventSet;

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

    double getMHAcceptanceRate();
    void resetMHAcceptanceParameters();

    void setCurrentLogLikelihood(double x);
    double getCurrentLogLikelihood();

    virtual double computeLogLikelihood() = 0;
    virtual double computeLogPrior() = 0;

    void initializeModelFromEventDataFile();

    int getNumberOfEvents();
    BranchEvent* getRootEvent();

    void changeNumberOfEventsMH();
    void moveEventMH();
    void updateEventRateMH();

    void addEventToTree();
    void addEventToTree(double x);

    void deleteEventFromTree(BranchEvent* be);
    void deleteRandomEventFromTree();

    void getEventDataString(std::stringstream& ss);

    double getModelTemperatureMH(void);
    void setModelTemperatureMH(double x);
    
protected:

    BranchEvent* chooseEventAtRandom();

    void addEventMH();
    void removeEventMH();

    void revertMovedEventToPrevious();
    void restoreLastDeletedEvent();

    double computeEventGainLogHR(double K, double logLikelihood,
        double oldLogLikelihood, double logPrior, double oldLogPrior,
            double qRatio);

    double computeEventLossLogHR(double K, double logLikelihood,
        double oldLogLikelihood, double logPrior, double oldLogPrior,
            double qRatio);

    void eventMove(bool local);
    void eventLocalMove();
    void eventGlobalMove();

    bool acceptMetropolisHastings(double lnR);
    bool isEventConfigurationValid(BranchEvent* be);

    // These functions take a branch event and recursively update
    // the branch histories for all nodes going toward the tips
    void forwardSetBranchHistories(BranchEvent* x);
    void forwardSetHistoriesRecursive(Node* p);

    double safeExponentiation(double x);

    // Pure virtual methods to be implemented by derived classes

    virtual void readModelSpecificParameters(std::ifstream& inputFile) = 0;
    virtual void setRootEventWithReadParameters() = 0;
    virtual BranchEvent* newBranchEventWithReadParameters
        (Node* x, double time) = 0;
    virtual void setMeanBranchParameters() = 0;

    virtual BranchEvent* newBranchEventWithRandomParameters(double x) = 0;

    virtual void setDeletedEventParameters(BranchEvent* be) = 0;
    virtual double calculateLogQRatioJump() = 0;

    virtual BranchEvent* newBranchEventFromLastDeletedEvent() = 0;

    virtual void getSpecificEventDataString
        (std::stringstream& ss, BranchEvent* event) = 0;

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

    double _logLikelihood;

    int _acceptCount;
    int _rejectCount;
    int _acceptLast;    // true if last generation was accept; false otherwise
    // 0 = last was rejected; 1 = accepted; -1 = not set.

    EventSet _eventCollection;
    BranchEvent* _rootEvent;

    double _lastDeletedEventMapTime;    // map time of last deleted event

    // Last event modified, whether it is moved, or has value updated
    BranchEvent* _lastEventModified;

    // This parameter holds the density of the new parameters proposed
    // during jump moves. If the parameters are sampled from the prior,
    // these should exactly cancel.
    double _logQRatioJump;

    // Temperature parameter for Metropolis coupling:
    double _temperatureMH;

    double computeLogHastingsRatio(double logLikRatio,
                                    double logPriorRatio, double logQratio);
    
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


inline int Model::getNumberOfEvents()
{
    return (int)_eventCollection.size();
}


inline BranchEvent* Model::getRootEvent()
{
    return _rootEvent;
}


inline void Model::setCurrentLogLikelihood(double x)
{
    _logLikelihood = x;
}


inline double Model::getCurrentLogLikelihood()
{
    return _logLikelihood;
}

inline double Model::getModelTemperatureMH()
{
    return _temperatureMH;
}


inline double Model::computeLogHastingsRatio(double logLikRatio,
                                            double logPriorRatio, double logQratio)
{
    return (_temperatureMH * (logLikRatio + logPriorRatio) + logQratio);
    
}




#endif
