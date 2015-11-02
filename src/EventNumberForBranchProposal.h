#ifndef EVENT_NUMBER_FOR_BRANCH_PROPOSAL_H
#define EVENT_NUMBER_FOR_BRANCH_PROPOSAL_H


#include "Proposal.h"

class Random;
class Settings;
class Model;
class BranchEvent;


class EventNumberForBranchProposal : public Proposal
{
    enum ProposalType {
        AddEvent,
        RemoveEvent
    };

public:

    EventNumberForBranchProposal
        (Random& random, Settings& settings, Model& model);

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

private:

    double computeLogLikelihoodRatio();
    double computeLogPriorRatio();
    double computeLogQRatio();

    Random& _random;
    Model& _model;

    bool _validateEventConfiguration;

    int _numberOfBranches;
    double _totalTreeLength;

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
