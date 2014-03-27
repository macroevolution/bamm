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

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

private:

    double computeLogLikelihoodRatio();
    double computeLogPriorRatio();
    double computeLogQRatio();

    MbRandom& _rng;
    Settings& _settings;
    Model& _model;
    Prior& _prior;

    double _updateEventRateScale;

    double _currentEventRate;
    double _proposedEventRate;

    double _cterm;
};


#endif
