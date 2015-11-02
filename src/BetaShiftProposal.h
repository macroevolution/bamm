#ifndef BETA_SHIFT_PROPOSAL_H
#define BETA_SHIFT_PROPOSAL_H


#include "EventParameterProposal.h"

class Random;
class Settings;
class Model;
class Prior;


class BetaShiftProposal : public EventParameterProposal
{
public:

    BetaShiftProposal(Random& random, Settings& settings, Model& model,
        Prior& prior);

    virtual double acceptanceRatio();

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
