#ifndef BETA_INIT_PROPOSAL_H
#define BETA_INIT_PROPOSAL_H


#include "EventParameterProposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;


class BetaInitProposal : public EventParameterProposal
{
public:

    BetaInitProposal(MbRandom& rng, Settings& settings, Model& model,
        Prior& prior);

private:

    virtual double getCurrentParameterValue();
    virtual double computeNewParameterValue();

    virtual void setProposedParameterValue();
    virtual void revertToOldParameterValue();

    virtual void updateParameterOnTree();

    virtual double computeRootLogPriorRatio();
    virtual double computeNonRootLogPriorRatio();
    virtual double computeLogQRatio();

    double _updateBetaInitScale;
    double _cterm;
};


#endif
