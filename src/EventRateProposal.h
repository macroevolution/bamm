#ifndef EVENT_RATE_PROPOSAL_H
#define EVENT_RATE_PROPOSAL_H


#include "Proposal.h"

class Random;
class Settings;
class Model;
class Prior;


class EventRateProposal : public Proposal
{
public:

    EventRateProposal(Random& random, Settings& settings, Model& model,
        Prior& prior);

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

private:

    double computeLogLikelihoodRatio();
    double computeLogPosteriorRatio();
    double computeLogPriorRatio();
    double computeLogQRatio();

    Random& _random;
    Settings& _settings;
    Model& _model;
    Prior& _prior;

    double _updateEventRateScale;

    double _currentEventRate;
    double _proposedEventRate;

    double _cterm;
};


#endif
