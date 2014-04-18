#ifndef MU_SHIFT_PROPOSAL_H
#define MU_SHIFT_PROPOSAL_H


#include "EventParameterProposal.h"

class Random;
class Settings;
class Model;
class Prior;


class MuShiftProposal : public EventParameterProposal
{
public:

    MuShiftProposal(Random& random, Settings& settings, Model& model,
        Prior& prior);

private:

    virtual double getCurrentParameterValue();
    virtual double computeNewParameterValue();

    virtual void setProposedParameterValue();
    virtual void revertToOldParameterValue();

    virtual void updateParameterOnTree();

    virtual double computeRootLogPriorRatio();
    virtual double computeNonRootLogPriorRatio();

    double _updateMuShiftScale;
};


#endif
