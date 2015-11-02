#ifndef BRANCH_HISTORY_H
#define BRANCH_HISTORY_H

#include <set>
#include "BranchEvent.h"


class BranchHistory
{
    typedef std::set<BranchEvent*, BranchEvent::PtrCompare> EventSet;
    typedef EventSet::size_type EventSetSizeType;

private:

    BranchEvent* _nodeEvent;             // event describing focal node
    BranchEvent* _ancestralNodeEvent;    // event describing ancestor

    // Set of all events on branch. Of length 0 if no events occurred on branch.
    // Also, if no events occur on branch, then entire branch is described by
    // the event referenced at nodeEvent
    EventSet _eventsOnBranch;

public:

    BranchHistory();

    BranchEvent* getLastEvent();

    // Get last event from a reference event that occurred on branch:
    BranchEvent* getLastEvent(BranchEvent* x);
    BranchEvent* getLastEvent(double ttime);
    BranchEvent* getNextEvent(double ttime);
    int getNumberOfEventsOnInterval(double t1, double t2);
    BranchEvent* getEventByIndexPosition(int i);

    void   setNodeEvent(BranchEvent* x);
    BranchEvent* getNodeEvent();

    void   setAncestralNodeEvent(BranchEvent* x);
    BranchEvent* getAncestralNodeEvent();

    void printBranchHistory();
    void reversePrintBranchHistory();
    void printEvent(BranchEvent* event);

    void popEventOffBranchHistory(BranchEvent* x);
    void addEventToBranchHistory(BranchEvent* x);
    int  getNumberOfBranchEvents();
};


#endif
