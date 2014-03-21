#ifndef LAMBDA_INIT_PROPOSAL_H
#define LAMBDA_INIT_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;
class SpExBranchEvent;


class LambdaInitProposal : public Proposal
{
public:

    LambdaInitProposal(MbRandom& rng, Settings& settings, Model& model,
        Prior& prior);

private:

    virtual void saveCurrentState();
    virtual void proposeNewState();

    double computeLogLikelihoodRatio();
    double computeLogPriorRatio();
    double computeLogQRatio();

    void specificAccept();
    void specificReject();

    Prior& _prior;

    SpExBranchEvent* _event;

    double _currentLambdaInit;
    double _proposedLambdaInit;

    double _currentLogLikelihood;
    double _proposedLogLikelihood;

    double _cterm;

    double _updateLambdaInitScale;
};


#endif
