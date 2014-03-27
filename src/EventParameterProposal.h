#ifndef EVENT_PARAMETER_PROPOSAL_H
#define EVENT_PARAMETER_PROPOSAL_H


#include "Proposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;
class Tree;
class BranchEvent;


class EventParameterProposal : public Proposal
{
public:

    EventParameterProposal(MbRandom& rng, Settings& settings, Model& model,
        Prior& prior);

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

protected:

    virtual double getCurrentParameterValue() = 0;
    virtual double computeNewParameterValue() = 0;

    virtual void setProposedParameterValue() = 0;
    virtual void revertToOldParameterValue() = 0;

    virtual double computeLogLikelihoodRatio();
    virtual double computeLogPriorRatio();
    virtual double computeRootLogPriorRatio();
    virtual double computeNonRootLogPriorRatio();
    virtual double computeLogQRatio();

    MbRandom& _rng;
    Settings& _settings;
    Model& _model;
    Prior& _prior;

    Tree* _tree;
    BranchEvent* _event;

    double _currentParameterValue;
    double _proposedParameterValue;

    double _currentLogLikelihood;
    double _proposedLogLikelihood;
};


#endif
