#ifndef BRANCH_EVENT_H
#define BRANCH_EVENT_H

#include "Tree.h"
#include "Log.h"

class Node;
class Random;


// This base class contains:
// (1) the node associated with the event
// (2) the map position of the event

class BranchEvent
{

public:

    class PtrCompare
    {
    public:
        bool operator()(const BranchEvent* m1,
                        const BranchEvent* m2) const {
            return *m1 < *m2;
        }
    };

private:

    double mapTime;
    Node* nodeptr;
    Tree* treePtr;
    Random& _random;

    // Keep values of the old pointer and old maptime associated
    // with the event for FAST reference if a LOCAL proposal is rejected.

    Node* oldNodePtr;
    double oldMapTime;

    double _absTime;    // real time, measured with t = 0 at root.

    // Allow rjMCMC to move between time-varying and time-constant partitions.
    bool _isEventTimeVariable;

public:

    BranchEvent(Node* x, Tree* tp, Random& random, double map);
    virtual ~BranchEvent();

    void   setMapTime(double x);
    double getMapTime();

    void  setEventNode(Node* x);
    Node* getEventNode();

    void   setAbsoluteTime(double x);
    double getAbsoluteTime();

    double getTimeToTip();

    void incrementMapPosition(double ink);
    void moveEventLocal(double stepsize);
    void moveEventGlobal();
    void setEventByMapPosition(double x);

    // Functions to set and manipulate old events:
    void  setOldEventNode(Node* x);
    Node* getOldEventNode();

    void   setOldMapTime(double x);
    double getOldMapTime();

    // Revert to old map position using oldPtr and oldMapTime
    // this only works if you have changed the nodeptr and maptime
    // relative to the values of oldNodePtr and oldMapTime
    void revertOldMapPosition();

    // Overloading comparision operator:
    bool operator<(const BranchEvent& a) const;

    // For time-varying rjMCMC:
    void setIsEventTimeVariable(bool x);
    bool getIsEventTimeVariable();
};


inline void BranchEvent::setMapTime(double x)
{
    mapTime = x;
}


inline double BranchEvent::getMapTime()
{
    return mapTime;
}


inline void BranchEvent::setEventNode(Node* x)
{
    nodeptr = x;
}


inline Node* BranchEvent::getEventNode()
{
    return nodeptr;
}


inline void BranchEvent::setAbsoluteTime(double x)
{
    _absTime = x;
}


inline double BranchEvent::getAbsoluteTime()
{
    return _absTime;
}


inline double BranchEvent::getTimeToTip()
{
    return treePtr->getTreeLength() - _absTime;
}


inline void BranchEvent::setOldEventNode(Node* x)
{
    oldNodePtr = x;
}


inline Node* BranchEvent::getOldEventNode()
{
    return oldNodePtr;
}


inline void BranchEvent::setOldMapTime(double x)
{
    oldMapTime = x;
}


inline double BranchEvent::getOldMapTime()
{
    return oldMapTime;
}


inline void BranchEvent::setIsEventTimeVariable(bool x)
{
    _isEventTimeVariable = x;
}


inline bool BranchEvent::getIsEventTimeVariable()
{
    return _isEventTimeVariable;
}


#endif
