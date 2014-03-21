#ifndef PROPOSAL_H
#define PROPOSAL_H


class MbRandom;
class Settings;
class Model;


class Proposal
{
public:

    Proposal(MbRandom& rng, Settings& settings, Model& model);
    virtual ~Proposal();

    void propose();
    void accept();
    void reject();

protected:

    virtual void saveCurrentState() = 0;
    virtual void proposeNewState() = 0;

    void updateLogRatios();
    virtual double computeLogLikelihoodRatio() = 0;
    virtual double computeLogPriorRatio() = 0;
    virtual double computeLogQRatio() = 0;

    virtual void specificAccept() = 0;
    virtual void specificReject() = 0;

    MbRandom& _rng;
    Settings& _settings;
    Model& _model;
};


#endif
