#ifndef SP_EX_MODEL_H
#define SP_EX_MODEL_H


#include "Model.h"
#include "LambdaInitProposal.h"
#include "LambdaShiftProposal.h"
#include "MuInitProposal.h"
#include "MuShiftProposal.h"

#include <iosfwd>
#include <vector>
#include <string>

class Node;
class Random;
class Settings;
class BranchEvent;
class Proposal;


class SpExModel : public Model
{

public:

    SpExModel(Random& rng, Settings& settings);

    virtual double computeLogLikelihood();
    virtual double computeLogPrior();

	// Methods for auto-tuning
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);

private:

    virtual void initializeSpecificUpdateWeights();

    virtual Proposal* getSpecificProposal(int parameter);
    
    virtual void setRootEventWithReadParameters
        (const std::vector<std::string>& parameters);
    virtual BranchEvent* newBranchEventWithReadParameters
        (Node* x, double time, const std::vector<std::string>& parameters);

    double lambdaInitParameter(const std::vector<std::string>& parameters);
    double lambdaShiftParameter(const std::vector<std::string>& parameters);
    double muInitParameter(const std::vector<std::string>& parameters);
    double muShiftParameter(const std::vector<std::string>& parameters);

    virtual BranchEvent* newBranchEventWithRandomParameters(double x);
    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    virtual void setMeanBranchParameters();
    virtual void setDeletedEventParameters(BranchEvent* be);

    double computeSpExProbBranch(Node* node);
    void computeSpExProb(double& spProb, double& exProb,
        double lambda, double mu, double D0, double E0, double deltaT);

    virtual double calculateLogQRatioJump();

    LambdaInitProposal _lambdaInitProposal;
    LambdaShiftProposal _lambdaShiftProposal;
    MuInitProposal _muInitProposal;
    MuShiftProposal _muShiftProposal;

    // Root event parameters
    double _lambdaInit0;
    double _lambdaShift0;
    double _muInit0;
    double _muShift0;

    bool _sampleFromPriorOnly;

    double _lastDeletedEventLambdaInit;
    double _lastDeletedEventLambdaShift;
    double _lastDeletedEventMuInit;
    double _lastDeletedEventMuShift;

    double _segLength;

    double _readLambdaInit;
    double _readLambdaShift;
    double _readMuInit;
    double _readMuShift;

    double _extinctionProbMax;
};


#endif
