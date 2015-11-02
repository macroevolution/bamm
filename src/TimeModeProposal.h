#ifndef TIME_MODE_PROPOSAL_H
#define TIME_MODE_PROPOSAL_H


#include "Proposal.h"
#include "Prior.h"

class Random;
class Settings;
class Model;
class Tree;
class BranchEvent;


class TimeModeProposal : public Proposal
{
    enum ProposalType {
        TimeConstant,
        TimeVariable
    };

public:

    TimeModeProposal(Random& random, Settings& settings, Model& model);

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

protected:

    virtual double initialParameter(BranchEvent* event) = 0;
    virtual double rateParameter(BranchEvent* event) = 0;
    virtual bool isTimeVariable(BranchEvent* event) = 0;
    
    virtual void setEventParameters(BranchEvent* event,
        double initParam, double rateParam, bool isTimeVariable) = 0;

    virtual void setModelParameters() = 0;

    virtual double rateParameterFromPrior() = 0;

    void makeTimeConstant(BranchEvent* event);
    void makeTimeVariable(BranchEvent* event);

    double computeMeanRate(double init, double k, double T);
    double computeRateInit(double mean, double k, double T);

    double computeLogLikelihoodRatio();
    double computeLogPriorRatio();
    double computeLogJacobian();
    double computeJacobian(double k, double T);

    Model& _model;
    Tree* _tree;

    Prior _prior;

    BranchEvent* _event;

    double _currentInitParam;
    double _currentRateParam;
    bool _currentIsTimeVariable;

    double _currentLogLikelihood;
    double _currentLogPrior;

    double _proposedLogLikelihood;
    double _proposedLogPrior;

    ProposalType _lastTimeModeProposal;
};


#endif
