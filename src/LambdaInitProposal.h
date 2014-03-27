#ifndef LAMBDA_INIT_PROPOSAL_H
#define LAMBDA_INIT_PROPOSAL_H


#include "EventParameterProposal.h"

class MbRandom;
class Settings;
class Model;
class Prior;


class LambdaInitProposal : public EventParameterProposal
{
public:

    LambdaInitProposal(MbRandom& rng, Settings& settings, Model& model,
        Prior& prior);

private:

    virtual double getCurrentParameterValue();
    virtual double computeNewParameterValue();

    virtual void setProposedParameterValue();
    virtual void revertToOldParameterValue();

    virtual double computeRootLogPriorRatio();
    virtual double computeNonRootLogPriorRatio();
    virtual double computeLogQRatio();

    double _updateLambdaInitScale;
    double _cterm;
};


#endif
