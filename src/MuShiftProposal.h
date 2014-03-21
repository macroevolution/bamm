#ifndef MU_SHIFT_PROPOSAL_H
#define MU_SHIFT_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;
class SpExBranchEvent;


class MuShiftProposal : public Proposal
{
public:

    MuShiftProposal(MbRandom& rng, Settings& settings, Model& model,
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

    double _currentMuShift;
    double _proposedMuShift;

    double _currentLogLikelihood;
    double _proposedLogLikelihood;

    double _updateMuShiftScale;
};


#endif
