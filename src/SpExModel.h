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
class SpExBranchEvent;


class SpExModel : public Model
{

public:

    SpExModel(Random& rng, Settings& settings);

    virtual double computeLogLikelihood();
    virtual double computeLogPrior();
 
	// Methods for auto-tuning
    //   no auto-tuning yet implemented
	void setUpdateLambdaInitScale(double x);
	void setUpdateMuInitScale(double x);
	void setUpdateLambdaShiftScale(double x);
    
    
    // FOSSIL STUFF
    double getPreservationRate();
    void setPreservationRate(double x);
    
    bool   getHasPaleoData();
    void   setHasPaleoData(bool x);
    
    void   initializeHasPaleoData();
    
    // Debugging likelihood function:
    void   printNodeProbs();

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
    virtual BranchEvent* newBranchEventWithParametersFromSettings(double x);
    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    virtual void setMeanBranchParameters();
    virtual void setDeletedEventParameters(BranchEvent* be);

    double computeSpExProbBranch(Node* node);
    void computeSpExProb(double& spProb, double& exProb,
        double lambda, double mu, double psi, double D0, double E0, double deltaT);

    
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
    int     _numberOccurrences;
    bool   _hasPaleoData;
    
    bool   _conditionOnSurvival;
    
    double _tmpvar;
    
    // October 2015 fixes
    bool _lambdaIsTimeVariable;
    
    
    // Debugging
    std::vector<double> _Dinitvec;
    std::vector<double> _Dfinalvec;
    std::vector<double> _Einitvec;
    std::vector<double> _Efinalvec;
    int nodeIndexLookup(Node* Node);
    
    void initializeDebugVectors();
    void outputDebugVectors();
    
     
    // BAMM updates October 2015:
    double computeMeanExponentialRateForInterval
                    (double rate_init, double rate_shift, double t_start, double t_end);
    double recomputeE0(double start_time, double end_time, double lam_init, double lam_shift,
                        double mu_init, double mu_shift, double Etip);
    
    bool _alwaysRecomputeE0;
    
    std::string _combineExtinctionAtNodes;
    
    
    
};


inline double SpExModel::getPreservationRate(void)
{
    return _preservationRate;
}

inline void SpExModel::setPreservationRate(double x)
{
    _preservationRate = x;
}


inline bool SpExModel::getHasPaleoData()
{
    return _hasPaleoData;
}


inline void SpExModel::setHasPaleoData(bool x)
{
    _hasPaleoData = x;
}



#endif







