#ifndef PROPOSAL_H
#define PROPOSAL_H


class Proposal
{
public:

    virtual ~Proposal();

    virtual void propose() = 0;
    virtual void accept() = 0;
    virtual void reject() = 0;

    virtual double acceptanceRatio() = 0;
};


#endif
