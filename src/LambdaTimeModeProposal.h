#ifndef LAMBDA_TIME_MODE_PROPOSAL_H
#define LAMBDA_TIME_MODE_PROPOSAL_H


#include "Proposal.h"
#include "Prior.h"

class Random;
class Settings;
class Model;
class Tree;
class SpExBranchEvent;


class LambdaTimeModeProposal : public Proposal
{
    enum ProposalType {
        TimeConstant,
        TimeVariable
    };

public:

    LambdaTimeModeProposal(Random& random, Settings& settings, Model& model);

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

private:

    void makeTimeConstant(SpExBranchEvent* event);
    void makeTimeVariable(SpExBranchEvent* event);

    double computeMeanLambda(double lam0, double k, double T);
    double computeLambdaInit(double lam, double k, double T);

    double computeLogLikelihoodRatio();
    double computeLogPriorRatio();
    double computeLogJacobian();
    double computeJacobian(double k, double T);

    Model& _model;
    Tree* _tree;

    Prior _prior;

    SpExBranchEvent* _event;

    double _currentLambdaInit;
    double _currentLambdaShift;
    bool _currentIsTimeVariable;

    double _currentLogLikelihood;
    double _currentLogPrior;

    double _proposedLogLikelihood;
    double _proposedLogPrior;

    ProposalType _lastTimeModeProposal;
};


#endif
