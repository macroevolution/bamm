#ifndef SP_EX_MODEL_H
#define SP_EX_MODEL_H

#include "Model.h"

#include <set>
#include <vector>
#include <sstream>

#include "BranchEvent.h"

//Forward declarations
class Tree;
class Node;
class MbRandom;
class Settings;
class Prior;


class SpExModel : public Model
{

public:

    SpExModel(MbRandom* rng, Tree* tree, Settings* settings, Prior* prior);
    virtual ~SpExModel();

    // Full likelihood will be lnLikTraits + lnLikBranches
    void   setCurrLnLTraits(double x);
    double getCurrLnLTraits();

    void   setCurrLnLBranches(double x);
    double getCurrLnLBranches();

    // Full likelihood function


    // Likelihood functions for branches data
    double computeLikelihoodBranches();
    double computeLikelihoodBranchesByInterval();

    double computeLogPrior();

    double safeExponentiation(double x);
    double proportionalShrink(double x, double scale);
    bool   acceptMetropolisHastings(const double lnR);

    // Stuff for event handling:
    void deleteRandomEventFromTree();
    void deleteEventFromTree(BranchEvent* be);

    int countEventsInBranchHistory(Node* p);

    // Initialize all branch histories to the root node.
    void initializeBranchHistories(Node* x);

    void printStartAndEndEventStatesForBranch(Node* x);

    void eventLocalMove(BranchEvent* x);   // move specific event
    void eventGlobalMove(BranchEvent* x);  // move specific event
    void eventLocalMove();             // move random event
    void eventGlobalMove();            // move random event
    void revertEventToPreviousPosition();

    // MCMC:

    // Propose addition or deletion; accept/reject move.
    void changeNumberOfEventsMH();
    void moveEventMH();
    void revertMovedEventToPrevious();

    /*  ***************** */
    // lambda/mu related stuff:
    void updateLambdaInitMH();
    void updateLambdaShiftMH();

    void updateMuInitMH();
    void updateMuShiftMH();

    void updateNodeStateMH();
    void updateNodeStateMH(Node* xnode);
    void updateDownstreamNodeStatesMH(Node* xnode);
    void updateEventRateMH();
    void updateTimeVariablePartitionsMH();

    /*  ********************  */
    void printEventData();

    void restoreLastDeletedEvent();

    // Troubleshooting
    void printBranchHistories(Node* x);

    // More output: acceptance rates
    double getMHacceptanceRate();
    void   resetMHacceptanceParameters();

    BranchEvent*  getEventByIndex(int x);

    void printExtinctionParams();
    int countTimeVaryingRatePartitions();

    // Generate string with event data:
    void getEventDataString(std::stringstream& ss);

    bool isEventConfigurationValid(BranchEvent* be);

    void debugLHcalculation();
	
	// Functions for auto-tuning
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);

private:

    virtual void readModelSpecificParameters(std::ifstream& inputFile);
    virtual void setRootEventWithReadParameters();
    virtual BranchEvent* newBranchEventWithReadParameters(Node* x, double time);
    virtual void setMeanBranchParameters();
    virtual BranchEvent* newBranchEventWithRandomParameters(double x);

    // Parameters of the model:

    double lnLikTraits;
    double lnLikBranches;

    // new parameters: March 23 2012
    double _updateLambdaInitScale;
    double _updateLambdaShiftScale;

    double _updateMuInitScale;
    double _updateMuShiftScale;

    // This parameter holds the density of the new
    // parameters proposed during jump moves.
    // If the parameters are sampled from the prior,
    //    these should exactly cancel.

    double _logQratioJump;

    // Root event parameters:

    double _lambdaInit0;
    double _lambdaShift0;

    double _muInit0;
    double _muShift0;

    // Priors
    double _lambdaInitPrior;
    double _lambdaShiftPrior;

    double _muInitPrior;
    double _muShiftPrior;

    double _lastDeletedEventLambdaInit;
    double _lastDeletedEventLambdaShift;

    double _lastDeletedEventMuInit;
    double _lastDeletedEventMuShift;

    double _segLength; // for splitting branches

    double _readLambdaInit;
    double _readLambdaShift;
    double _readMuInit;
    double _readMuShift;
};


inline void SpExModel::setCurrLnLTraits(double x)
{
    lnLikTraits = x;
}


inline double SpExModel::getCurrLnLTraits()
{
    return lnLikTraits;
}


inline void SpExModel::setCurrLnLBranches(double x)
{
    lnLikBranches = x;
}


inline double SpExModel::getCurrLnLBranches()
{
    return lnLikBranches;
}


inline void SpExModel::setUpdateLambdaInitScale(double x)
{
	_updateLambdaInitScale = x;
}


inline void SpExModel::setUpdateMuInitScale(double x)
{
	_updateMuInitScale = x;
}


inline void SpExModel::setUpdateLambdaShiftScale(double x)
{
	_updateLambdaShiftScale = x;
}


#endif
