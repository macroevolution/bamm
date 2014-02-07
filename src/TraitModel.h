#ifndef TRAIT_MODEL_H
#define TRAIT_MODEL_H

#include "Model.h"

#include <set>
#include <vector>
#include <sstream>

#include "BranchEvent.h"

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
    void addEventToTreeWithSetBeta(double beta, double bshift);
    void deleteRandomEventFromTree();
    void deleteEventFromTree(BranchEvent* be);

    int countEventsInBranchHistory(Node* p);

    // initialize all branch histories to the root node.
    void initializeBranchHistories(Node* x);
    void printStartAndEndEventStatesForBranch(Node* x);

    // MCMC:
    // Propose addition or deletion; accept/reject move.
    void changeNumberOfEventsMH();
    void moveEventMH();

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

    BranchEvent* getEventByIndex(int x);

    int countTimeVaryingRatePartitions();

    // Generate string with event data:
    void getEventDataString(std::stringstream& ss);
    bool isEventConfigurationValid(BranchEvent* be);

    double      getLastLH();

private:

    virtual void readModelSpecificParameters(std::ifstream& inputFile);
    virtual void setRootEventWithReadParameters();
    virtual BranchEvent* newBranchEventWithReadParameters(Node* x, double time);
    virtual void setMeanBranchParameters();

    virtual BranchEvent* newBranchEventWithRandomParameters(double x);

    //  parameters of the model:

    double lnLikTraits;

    double _updateBetaScale;
    double _updateBetaShiftScale;
    double _updateNodeStateScale;

    double _lastDeletedEventBetaInit;;
    double _lastDeletedEventBetaShift;

    // Here are several variables that track the previous
    // state. At some point, these should have their own class

    // Ultimately, initializations should be handled in TraitModel
    // and Model classes
    // NOT in class TREE!!
    void   initializeTraitParamsForNodes();
    double _lastLH;

    double _logQratioJump;

    double _readBetaInit;
    double _readBetaShift;
};


inline void TraitModel::setCurrLnLTraits(double x)
{
    lnLikTraits = x;
}


inline double TraitModel::getCurrLnLTraits()
{
    return lnLikTraits;
}


inline double TraitModel::getLastLH()
{
    return _lastLH;
}


#endif
