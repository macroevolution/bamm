#ifndef BETA_SHIFT_PROPOSAL_H
#define BETA_SHIFT_PROPOSAL_H


#include "EventParameterProposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;


class BetaShiftProposal : public EventParameterProposal
{
public:

    BetaShiftProposal(MbRandom& rng, Settings& settings, Model& model,
        Prior& prior);

private:

    virtual double getCurrentParameterValue();
    virtual double computeNewParameterValue();

    virtual void setProposedParameterValue();
    virtual void revertToOldParameterValue();

    virtual void updateParameterOnTree();

    virtual double computeRootLogPriorRatio();
    virtual double computeNonRootLogPriorRatio();

    double _updateBetaShiftScale;
};


#endif
