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

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

private:

    virtual double computeLogLikelihoodRatio();

    MbRandom& _rng;
    Settings& _settings;
    Model& _model;

    double _localToGlobalMoveRatio;
    double _scale;

    bool _validateEventConfiguration;

    BranchEvent* _event;

    int _currentEventCount;
    double _currentLogLikelihood;
    double _proposedLogLikelihood;
};


#endif
