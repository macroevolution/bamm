/*
 *  BranchHistory.h
 *  rateshift
 *
 *  Created by Dan Rabosky on 12/14/11.
 *
 ptr to BranchHistory will be included with Node constructor
 each node has reference to all events associated with the branch
 If no events occur on branch, you still can obtain the event reference for
 the node itself and the end of the branch very easily

 */


#ifndef BRANCHHISTORY_H
#define BRANCHHISTORY_H

#include <set>
#include "BranchEvent.h"


class BranchHistory
{

    typedef std::set<BranchEvent*, BranchEventPtrCompare> EventSet;
    typedef EventSet::size_type EventSetSizeType;

private:

    BranchEvent* nodeEvent;          // event describing focal node
    BranchEvent* ancestralNodeEvent; // event describing ancestor

    // Set of all events on branch. Of length 0 if no events occurred on branch.
    // Also, if no events occur on branch, then entire branch is described by
    // the event referenced at nodeEvent;
    EventSet eventsOnBranch;

public:

    BranchHistory();
    ~BranchHistory();

    BranchEvent* getLastEvent(); // get last event on branch

    // Get last event from a reference event that occurred on branch:
    BranchEvent* getLastEvent(BranchEvent* x);
    BranchEvent* getLastEvent(double ttime);
    BranchEvent* getNextEvent(double ttime);
    int getNumberOfEventsOnInterval(double t1, double t2);
    BranchEvent* getEventByIndexPosition(int i);

    void         setNodeEvent(BranchEvent* x);
    BranchEvent* getNodeEvent();

    void         setAncestralNodeEvent(BranchEvent* x);
    BranchEvent* getAncestralNodeEvent();

    void printBranchHistory();
    void reversePrintBranchHistory();

    void popEventOffBranchHistory(BranchEvent* x);
    void addEventToBranchHistory(BranchEvent* x);
    int  getNumberOfBranchEvents();
};


inline void BranchHistory::setNodeEvent(BranchEvent* x)
{
    nodeEvent = x;
}


inline BranchEvent* BranchHistory::getNodeEvent()
{
    return nodeEvent;
}


inline void BranchHistory::setAncestralNodeEvent(BranchEvent* x)
{
    ancestralNodeEvent = x;
}


inline BranchEvent* BranchHistory::getAncestralNodeEvent()
{
    return ancestralNodeEvent;
}


inline void BranchHistory::popEventOffBranchHistory(BranchEvent* x)
{
    eventsOnBranch.erase(x);
}


inline void BranchHistory::addEventToBranchHistory(BranchEvent* x)
{
    eventsOnBranch.insert(x);
}


inline int BranchHistory::getNumberOfBranchEvents()
{
    return (int)eventsOnBranch.size();
}


#endif
