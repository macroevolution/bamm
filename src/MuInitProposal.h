#ifndef MU_INIT_PROPOSAL_H
#define MU_INIT_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;
class SpExBranchEvent;


class MuInitProposal : public Proposal
{
public:

    MuInitProposal(MbRandom& rng, Settings& settings, Model& model,
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

    double _currentMuInit;
    double _proposedMuInit;

    double _currentLogLikelihood;
    double _proposedLogLikelihood;

    double _cterm;

    double _updateMuInitScale;
};


#endif
