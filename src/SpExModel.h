#ifndef SP_EX_MODEL_H
#define SP_EX_MODEL_H

#include "Model.h"
#include "LambdaInitProposal.h"
#include "LambdaShiftProposal.h"
#include "MuInitProposal.h"
#include "MuShiftProposal.h"
#include "LambdaTimeModeProposal.h"

#include <iosfwd>

class Node;
class Random;
class Settings;
class BranchEvent;
class Proposal;


class SpExModel : public Model
{

public:

    SpExModel(Random& rng, Settings* settings);

    virtual double computeLogLikelihood();
    virtual double computeLogPrior();

	// Methods for auto-tuning
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);

private:

    double computeLogLikelihoodByInterval();

    virtual void initializeSpecificUpdateWeights();

    virtual Proposal* getSpecificProposal(int parameter);
    
    virtual void readModelSpecificParameters(std::ifstream& inputFile);
    virtual void setRootEventWithReadParameters();

    virtual BranchEvent* newBranchEventWithReadParameters(Node* x, double time);
    virtual BranchEvent* newBranchEventWithRandomParameters(double x);
    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    virtual void setMeanBranchParameters();
    virtual void setDeletedEventParameters(BranchEvent* be);

    virtual double calculateLogQRatioJump();

    LambdaInitProposal _lambdaInitProposal;
    LambdaShiftProposal _lambdaShiftProposal;
    MuInitProposal _muInitProposal;
    MuShiftProposal _muShiftProposal;
    LambdaTimeModeProposal _lambdaTimeModeProposal;

    // Root event parameters
    double _lambdaInit0;
    double _lambdaShift0;
    double _muInit0;
    double _muShift0;
    bool _initialLambdaIsTimeVariable;

    bool _sampleFromPriorOnly;

    double _lastDeletedEventLambdaInit;
    double _lastDeletedEventLambdaShift;
    double _lastDeletedEventMuInit;
    double _lastDeletedEventMuShift;
    double _lastDeletedEventTimeVariable;

    double _segLength;

    double _readLambdaInit;
    double _readLambdaShift;
    double _readMuInit;
    double _readMuShift;

    double _extinctionProbMax;
};


#endif
