#ifndef SP_EX_MODEL_H
#define SP_EX_MODEL_H


#include "Model.h"

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

    
    double computePreservationLogProb();
    
    
    virtual double calculateLogQRatioJump();

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
    
    //FOSSIL
    // Fossil preservation rate. Assume 1 value for now.
    double _preservationRate;
    // Time at which tree is observed, relative to root
    double _observationTime;
    
    // this will change at some point
    //   to allow greater flexibility of preservation model
    double _numberOccurrences;
    
};


#endif







