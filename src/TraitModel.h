#ifndef TRAIT_MODEL_H
#define TRAIT_MODEL_H

#include "Model.h"

#include <set>
#include <vector>
#include <sstream>

#include "TraitBranchEvent.h"

// Forward declarations
class Tree;
class Node;
class MbRandom;
class Settings;
class Prior;


class TraitModel : public Model
{

public:

    TraitModel(MbRandom* rng, Tree* tree, Settings* settings, Prior* prior);
    virtual ~TraitModel();

    // Full likelihood will be lnLikTraits + lnLikBranches
    void   setCurrLnLTraits(double x);
    double getCurrLnLTraits();

    double computeLikelihoodTraits();
    double computeTriadLikelihoodTraits(Node* x);

    double computeLogPrior();

    double safeExponentiation(double x);
    double proportionalShrink(double x, double scale);
    bool   acceptMetropolisHastings(const double lnR);

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

    // more output: acceptance rates
    double getMHacceptanceRate();
    void   resetMHacceptanceParameters();

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

    double _updateBetaScale;
    double _updateBetaShiftScale;
    double _updateNodeStateScale;

    // Other private variables

    // Event collection does not contain the root event
    std::set<TraitBranchEvent*, TraitBranchEventPtrCompare> eventCollection;
    TraitBranchEvent* rootEvent; //branch event at root node; can't be modified

    double _lastDeletedEventBetaInit;;
    double _lastDeletedEventBetaShift;

    // Here are several variables that track the previous
    // state. At some point, these should have their own class

    // this is a pointer to the last event modified, whether
    // it is moved, or has lambda updated, or whatever.
    TraitBranchEvent* lastEventModified;

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


inline int TraitModel::getNumberOfEvents()
{
    return (int)eventCollection.size();
}


inline TraitBranchEvent* TraitModel::getRootEvent()
{
    return rootEvent;
}


inline double TraitModel::getLastLH()
{
    return _lastLH;
}


#endif
