/*
 *  TraitModel.h
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.

 *
 */

#ifndef TRAITMODEL_H
#define TRAITMODEL_H

#include <set>
#include <vector>
#include <sstream>

//Forward declarations
class Tree;
class Node;
class MbRandom;
class Settings;
class TraitBranchEvent;

double safeExponentiation(double x);
double proportionalShrink(double x, double scale);


class TraitModel
{

public:

    static double mhColdness;

    TraitModel(MbRandom* ranptr, Tree* tp, Settings* sp);
    ~TraitModel();

    // Full likelihood will be lnLikTraits + lnLikBranches
    void   setCurrLnLTraits(double x);
    double getCurrLnLTraits();

    double computeLikelihoodTraits();
    double computeTriadLikelihoodTraits(Node* x);

    double computeLogPrior();

    double getEventRate();
    void   setEventRate(double x);

    double proportionalShrink(double x, double scale);
    bool   acceptMetropolisHastings(const double lnR);
    void   incrementGeneration();
    int    getGeneration();
    void   resetGeneration();

    // Stuff for event handling:
    void addEventToTree();
    void addEventToTreeWithSetBeta(double beta, double bshift);
    void addEventToTree(double x); // add to specific map position...
    void deleteRandomEventFromTree();
    void deleteEventFromTree(TraitBranchEvent* be);

    void printEvents();

    int getNumberOfEvents();

    TraitBranchEvent* getRootEvent();

    // These functions take a branch event
    //  and recursively update branch histories for all nodes
    //  going towards the tips
    void forwardSetBranchHistories(TraitBranchEvent* x);
    void forwardSetHistoriesRecursive(Node* p);

    int countEventsInBranchHistory(Node* p);

    // initialize all branch histories to the root node.
    void initializeBranchHistories(Node* x);
    void printStartAndEndEventStatesForBranch(Node* x);


    void eventLocalMove(TraitBranchEvent* x); // move specific event
    void eventGlobalMove(TraitBranchEvent* x); // move specific event
    void eventLocalMove(); // move random event
    void eventGlobalMove(); // move random event
    void revertEventToPreviousPosition();

    // Return random event, or NULL if no events on tree (other than root)
    TraitBranchEvent* chooseEventAtRandom();

    // MCMC:
    // Propose addition or deletion; accept/reject move.
    void changeNumberOfEventsMH();
    void moveEventMH();
    void revertMovedEventToPrevious();

    /*  ***************** */

    // Trait evolution stuff:
    void updateBetaMH();
    void updateNodeStateMH();
    void updateNodeStateMH(Node* xnode);
    void updateBetaShiftMH();
    void updateDownstreamNodeStatesMH(Node* xnode);
    void updateEventRateMH();
    void updateTimeVariablePartitionsMH();
    void setMinMaxTraitPriors();

    /*  ********************  */

    void printEventData();
    void restoreLastDeletedEvent();

    // Troubleshooting
    void printBranchHistories(Node* x);

    Tree* getTreePtr();

    // more output: acceptance rates
    double getMHacceptanceRate();

    void setAcceptLastUpdate(int x);
    int  getAcceptLastUpdate();

    // 0 = last was rejected; 1 = accepted; -1 = not set.

    void   setPoissonRatePrior(double x);
    double getPoissonRatePrior();

    TraitBranchEvent* getEventByIndex(int x);

    int countTimeVaryingRatePartitions();

    // Generate string with event data:
    void getEventDataString(std::stringstream& ss);
    bool isEventConfigurationValid(TraitBranchEvent* be);

    double      getLastLH();

    void      initializeModelFromEventDataFileTrait();

private:

    //  parameters of the model:

    double lnLikTraits;

    double eventLambda; // Poisson rate

    double poissonRatePrior; // exponential  /* keep this in Model */

    double _updateBetaScale;
    double _updateBetaShiftScale;
    double _updateNodeStateScale;
    double _scale;
    double _updateEventRateScale;
    double _localGlobalMoveRatio;
    double _targetNumber;

    // Other private variables

    int gen;
    MbRandom* ran;
    Tree* treePtr;
    Settings* sttings;

    int acceptCount;
    int rejectCount;

    std::set<TraitBranchEvent*> eventCollection; // does NOT contain root event
    TraitBranchEvent* rootEvent; //branch event at root node; can't be modified

    double lastDeletedEventMapTime; // map time of last deleted event

    double _lastDeletedEventBetaInit;;
    double _lastDeletedEventBetaShift;

    // Here are several variables that track the previous
    // state. At some point, these should have their own class

    // this is a pointer to the last event modified, whether
    // it is moved, or has lambda updated, or whatever.
    TraitBranchEvent* lastEventModified;

    // General acceptreject flag:
    int acceptLast; // true if last generation was accept; false otherwise

    // Ultimately, initializations should be handled in TraitModel
    // and Model classes
    // NOT in class TREE!!
    void   initializeTraitParamsForNodes();
    double _lastLH;

    double _logQratioJump;
};


inline void TraitModel::setCurrLnLTraits(double x)
{
    lnLikTraits = x;
}


inline double TraitModel::getCurrLnLTraits()
{
    return lnLikTraits;
}


inline double TraitModel::getEventRate()
{
    return eventLambda;
}


inline void TraitModel::setEventRate(double x)
{
    eventLambda = x;
}


inline void TraitModel::incrementGeneration()
{
    gen++;
}


inline int TraitModel::getGeneration()
{
    return gen;
}


inline void TraitModel::resetGeneration()
{
    gen = 0;   // to be used after TraitPreBurnin
}


inline int TraitModel::getNumberOfEvents()
{
    return eventCollection.size();
}


inline TraitBranchEvent* TraitModel::getRootEvent()
{
    return rootEvent;
}


inline Tree* TraitModel::getTreePtr()
{
    return treePtr;
}


inline int TraitModel::getAcceptLastUpdate()
{
    return acceptLast;
}


inline void TraitModel::setAcceptLastUpdate(int x)
{
    acceptLast = x;
}


inline void TraitModel::setPoissonRatePrior(double x)
{
    poissonRatePrior  = x;
}


inline double TraitModel::getPoissonRatePrior()
{
    return poissonRatePrior;
}


inline double TraitModel::getLastLH()
{
    return _lastLH;
}


#endif
