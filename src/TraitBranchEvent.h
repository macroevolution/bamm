#ifndef TRAIT_BRANCH_EVENT_H
#define TRAIT_BRANCH_EVENT_H

#include "BranchEvent.h"

class Tree;
class Node;
class MbRandom;


class TraitBranchEvent : public BranchEvent
{

private:

    double _betaInit;  // initial beta value.
    double _betaShift; // temporal shift parameter of trait evolution rate.

public:

    // constructors, depending on whether you want trait rate or lambda/mu
    TraitBranchEvent(double beta, double shift,
            Node* x, Tree* tp, MbRandom* rp, double map);
    virtual ~TraitBranchEvent() {};

    void   setBetaInit(double x);
    double getBetaInit();

    void   setBetaShift(double x);
    double getBetaShift();

};


inline void TraitBranchEvent::setBetaInit(double x)
{
    _betaInit = x;
}


inline double TraitBranchEvent::getBetaInit()
{
    return _betaInit;
}


inline void TraitBranchEvent::setBetaShift(double x)
{
    _betaShift = x;
}


inline double TraitBranchEvent::getBetaShift()
{
    return _betaShift;
}


#endif
