#ifndef MOVE_EVENT_PROPOSAL_H
#define MOVE_EVENT_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Settings;
class Model;
class BranchEvent;


class MoveEventProposal : public Proposal
{
public:

    MoveEventProposal(MbRandom& rng, Settings& settings, Model& model);

private:

    virtual void saveCurrentState();
    virtual void proposeNewState();

    virtual double computeLogLikelihoodRatio();
    virtual double computeLogPriorRatio();
    virtual double computeLogQRatio();

    virtual void specificAccept();
    virtual void specificReject();

    double _localToGlobalMoveRatio;
    double _scale;

    BranchEvent* _event;

    int _currentEventCount;
    double _currentLogLikelihood;
    double _proposedLogLikelihood;
};


#endif
