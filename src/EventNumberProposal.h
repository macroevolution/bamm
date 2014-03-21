#ifndef EVENT_NUMBER_PROPOSAL_H
#define EVENT_NUMBER_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Model;
class BranchEvent;


class EventNumberProposal : public Proposal
{
    // TODO: Future: Add and Remove events should be separate Proposals
    enum ProposalType {
        AddEvent,
        RemoveEvent
    };

public:

    EventNumberProposal(MbRandom& rng, Settings& settings, Model& model);

private:

    virtual void saveCurrentState();
    virtual void proposeNewState();
    void proposeAddEvent();
    void proposeRemoveEvent();

    double computeLogLikelihoodRatio();
    double computeLogPriorRatio();
    double computeLogQRatio();

    void specificAccept();
    void specificReject();
    void rejectAddEvent();
    void rejectRemoveEvent();

    int _currentEventCount;
    double _currentLogLikelihood;
    double _currentLogPrior;

    int _proposedEventCount;
    double _proposedLogLikelihood;
    double _proposedLogPrior;

    ProposalType _lastProposal;
    BranchEvent* _lastEventChanged;
};


#endif
