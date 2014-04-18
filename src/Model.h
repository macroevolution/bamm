#ifndef MODEL_H
#define MODEL_H


#include "Prior.h"
#include "EventNumberProposal.h"
#include "MoveEventProposal.h"
#include "EventRateProposal.h"
#include "BranchEvent.h"

#include <vector>
#include <set>
#include <iosfwd>

class Random;
class Tree;
class Settings;
class Node;
class Proposal;


typedef std::set<BranchEvent*, BranchEvent::PtrCompare> EventSet;

class Model
{

public:

    Model(Random& random, Settings* settings);
    virtual ~Model();

    void finishConstruction();

    Tree* getTreePtr();

    double getEventRate();
    void setEventRate(double x);

    int getLastParameterUpdated();

    int getAcceptLastUpdate();
    void setAcceptLastUpdate(int x);

    double getMHAcceptanceRate();
    void resetMHAcceptanceParameters();

    void setCurrentLogLikelihood(double x);
    double getCurrentLogLikelihood();

    void setProposedLogLikelihood(double proposedLogLikelihood);

    virtual double computeLogLikelihood() = 0;
    virtual double computeLogPrior() = 0;

    void setLogLikelihoodRatio(double logLikelihoodRatio);
    void setLogPriorRatio(double logPriorRatio);
    void setLogQRatio(double logQRatio);

    void initializeModelFromEventDataFile();

    int getNumberOfEvents();
    BranchEvent* getRootEvent();

    EventSet& events();

    void proposeNewState();
    void acceptProposal();
    void rejectProposal();

    BranchEvent* chooseEventAtRandom(bool includeRoot = false);

    // These functions take a branch event and recursively update
    // the branch histories for all nodes going toward the tips
    void forwardSetBranchHistories(BranchEvent* x);
    void forwardSetHistoriesRecursive(Node* p);

    BranchEvent* addRandomEventToTree();
    BranchEvent* addEventToTree(BranchEvent* newEvent);

    BranchEvent* removeEventFromTree(BranchEvent* be);
    BranchEvent* removeRandomEventFromTree();

    virtual void setMeanBranchParameters() = 0;

    double getTemperatureMH(void);
    void setTemperatureMH(double x);

    double logQRatioJump();

    double acceptanceRatio();

    bool isEventConfigurationValid(BranchEvent* be);
    
protected:

    void calculateUpdateWeights();
    void initializeUpdateWeights();
    virtual void initializeSpecificUpdateWeights() = 0;

    int chooseParameterToUpdate();

    virtual Proposal* getSpecificProposal(int parameter) = 0;

    double computeLogHastingsRatio(double logLikRatio,
        double logPriorRatio, double logQratio);

    bool acceptMetropolisHastings(double lnR);

    double safeExponentiation(double x);

    // Pure virtual methods to be implemented by derived classes

    virtual void readModelSpecificParameters(std::ifstream& inputFile) = 0;
    virtual void setRootEventWithReadParameters() = 0;
    virtual BranchEvent* newBranchEventWithReadParameters
        (Node* x, double time) = 0;

    virtual BranchEvent* newBranchEventWithRandomParameters(double x) = 0;

    virtual void setDeletedEventParameters(BranchEvent* be) = 0;
    virtual double calculateLogQRatioJump() = 0;

    virtual BranchEvent* newBranchEventFromLastDeletedEvent() = 0;

    Random& _random;
    Settings* _settings;

    Prior _prior;

    Tree* _tree;

    EventNumberProposal _eventNumberProposal;
    MoveEventProposal _moveEventProposal;
    EventRateProposal _eventRateProposal;

    std::vector<double> _updateWeights;
    int _lastParameterUpdated;

    Proposal* _lastProposal;

    double _eventRate;    // Poisson rate

    double _logLikelihood;

    double _logLikelihoodRatio;
    double _logPriorRatio;
    double _logQRatio;

    double _proposedLogLikelihood;

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
};


inline Tree* Model::getTreePtr()
{
    return _tree;
}


inline int Model::getLastParameterUpdated()
{
    return _lastParameterUpdated;
}


inline EventSet& Model::events()
{
    return _eventCollection;
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


inline void Model::setLogLikelihoodRatio(double logLikelihoodRatio)
{
    _logLikelihoodRatio = logLikelihoodRatio;
}


inline void Model::setLogPriorRatio(double logPriorRatio)
{
    _logPriorRatio = logPriorRatio;
}


inline void Model::setLogQRatio(double logQRatio)
{
    _logQRatio = logQRatio;
}


inline void Model::setProposedLogLikelihood(double proposedLogLikelihood)
{
    _proposedLogLikelihood = proposedLogLikelihood;
}


inline double Model::getTemperatureMH()
{
    return _temperatureMH;
}


inline double Model::logQRatioJump()
{
    return _logQRatioJump;
}


#endif
