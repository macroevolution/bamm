#ifndef MU_SHIFT_PROPOSAL_H
#define MU_SHIFT_PROPOSAL_H


#include "EventParameterProposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;


class MuShiftProposal : public EventParameterProposal
{
public:

    MuShiftProposal(MbRandom& rng, Settings& settings, Model& model,
        Prior& prior);

private:

    virtual double getCurrentParameterValue();
    virtual double computeNewParameterValue();

    virtual void setProposedParameterValue();
    virtual void revertToOldParameterValue();

    virtual double computeRootLogPriorRatio();
    virtual double computeNonRootLogPriorRatio();

    double _updateMuShiftScale;
};


#endif
