#ifndef EVENT_RATE_PROPOSAL_H
#define EVENT_RATE_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;


class EventRateProposal : public Proposal
{
public:

    EventRateProposal(MbRandom& rng, Settings& settings, Model& model,
        Prior& prior);

private:

    virtual void saveCurrentState();
    virtual void proposeNewState();

    virtual double computeLogLikelihoodRatio();
    virtual double computeLogPriorRatio();
    virtual double computeLogQRatio();

    void specificAccept();
    void specificReject();

    Prior& _prior;

    double _updateEventRateScale;
    double _currentEventRate;
    double _proposedEventRate;

    double _cterm;
};


#endif
