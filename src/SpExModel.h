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

    // Likelihood functions for branches data
    virtual double computeLogLikelihood();
    virtual double computeLogPrior();

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

    void printExtinctionParams();
    int countTimeVaryingRatePartitions();

    // Generate string with event data:
    void getEventDataString(std::stringstream& ss);

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

    virtual void setDeletedEventParameters(BranchEvent* be);
    virtual double calculateLogQRatioJump();

    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    double computeLogLikelihoodByInterval();

    double _updateLambdaInitScale;
    double _updateLambdaShiftScale;

    double _updateMuInitScale;
    double _updateMuShiftScale;

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
