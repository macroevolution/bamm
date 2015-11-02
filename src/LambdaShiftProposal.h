#ifndef LAMBDA_SHIFT_PROPOSAL_H
#define LAMBDA_SHIFT_PROPOSAL_H


#include "EventParameterProposal.h"

class Random;
class Settings;
class Model;
class Prior;


class LambdaShiftProposal : public EventParameterProposal
{
public:

    LambdaShiftProposal(Random& random, Settings& settings, Model& model,
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

    double _updateLambdaShiftScale;
};


#endif
