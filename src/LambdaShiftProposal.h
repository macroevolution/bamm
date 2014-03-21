#ifndef LAMBDA_SHIFT_PROPOSAL_H
#define LAMBDA_SHIFT_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;
class SpExBranchEvent;


class LambdaShiftProposal : public Proposal
{
public:

    LambdaShiftProposal(MbRandom& rng, Settings& settings, Model& model,
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

    double _currentLambdaShift;
    double _proposedLambdaShift;

    double _currentLogLikelihood;
    double _proposedLogLikelihood;

    double _updateLambdaShiftScale;
};


#endif
