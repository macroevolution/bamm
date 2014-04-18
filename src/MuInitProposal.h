#ifndef MU_INIT_PROPOSAL_H
#define MU_INIT_PROPOSAL_H


#include "EventParameterProposal.h"

class Random;
class Settings;
class Model;
class Prior;


class MuInitProposal : public EventParameterProposal
{
public:

    MuInitProposal(Random& random, Settings& settings, Model& model,
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

    double _updateMuInitScale;
    double _cterm;
};


#endif
