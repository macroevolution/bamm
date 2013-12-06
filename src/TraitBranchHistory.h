#ifndef TRAITBRANCHHISTORY_H
#define TRAITBRANCHHISTORY_H

#include <set>
#include "TraitBranchEvent.h"


class TraitBranchHistory
{

    typedef std::set<TraitBranchEvent*, TraitBranchEventPtrCompare> EventSet;
    typedef EventSet::size_type EventSetSizeType;

private:

    TraitBranchEvent* nodeEvent; // event describing focal node
    TraitBranchEvent* ancestralNodeEvent; // event describing ancestor

    // set of all events on branch. Of length 0 if no events occurred on branch.
    // Also, if no events occur on branch, then entire branch is described by
    // the event referenced at nodeEvent;
    EventSet eventsOnBranch;

public:

    TraitBranchHistory();
    ~TraitBranchHistory();

    TraitBranchEvent* getLastEvent(); // get last event on branch

    // get last event from a reference event that occurred on branch:
    TraitBranchEvent* getLastEvent(TraitBranchEvent* x);
    TraitBranchEvent* getEventByIndexPosition(int i);

    void              setNodeEvent(TraitBranchEvent* x);
    TraitBranchEvent* getNodeEvent();

    void              setAncestralNodeEvent(TraitBranchEvent* x);
    TraitBranchEvent* getAncestralNodeEvent();

    void printBranchHistory();
    void reversePrintBranchHistory();

    void popEventOffBranchHistory(TraitBranchEvent* x);
    void addEventToBranchHistory(TraitBranchEvent* x);
    int  getNumberOfBranchEvents();
};


inline void TraitBranchHistory::setNodeEvent(TraitBranchEvent* x)
{
    nodeEvent = x;
}


inline TraitBranchEvent* TraitBranchHistory::getNodeEvent()
{
    return nodeEvent;
}


inline void TraitBranchHistory::setAncestralNodeEvent(TraitBranchEvent* x)
{
    ancestralNodeEvent = x;
}


inline TraitBranchEvent* TraitBranchHistory::getAncestralNodeEvent()
{
    return ancestralNodeEvent;
}


inline void TraitBranchHistory::popEventOffBranchHistory(TraitBranchEvent* x)
{
    eventsOnBranch.erase(x);
}


inline void TraitBranchHistory::addEventToBranchHistory(TraitBranchEvent* x)
{
    eventsOnBranch.insert(x);
}


inline int TraitBranchHistory::getNumberOfBranchEvents()
{
    return (int)eventsOnBranch.size();
}


#endif
