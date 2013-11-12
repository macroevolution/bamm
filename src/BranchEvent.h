/*
 *  BranchEvent.h
 *  rateshift
 *
 *  Created by Dan Rabosky on 12/5/11.


 // Major update no more phenotypic evolution
 // March 25 2012

 */

#ifndef BRANCHEVENT_H
#define BRANCHEVENT_H


//Forward declarations:
class Tree;
class Node;
class MbRandom;
class TraitBranchEvent;


/*

 class BranchEvent is the actual event.
 It contains:
  (1) the node associated with the event
  (2) the map position of the event
  (3) parameters associated with the event
    e.g., lambda
The initial event will always be the root node
 with a map position of 0.
This event is immutable.

 */

class BranchEvent
{

private:

    double mapTime;
    Node* nodeptr;
    Tree* treePtr;
    MbRandom* ranPtr;
    
    double epsilon;    // epsilon for local move

    // Keep values of the old pointer and old maptime associated
    // with the event for FAST reference if a LOCAL proposal is rejected.

    Node* oldNodePtr;
    double oldMapTime;

    /******************/
    // New event parameters, March 25, 2012
    double _absTime;   // real time, measured with t = 0 at root.
    double _lamInit;   // Initial speciation rate at event
    double _muInit;    // Initial Mu rate at event
    double _lamShift;  // magnitude & direction of speciation shift
    double _muShift;   // magnitude & direction of mu shift

    /*****************/
    // New parameters June 12 2012
    // allow rjMCMC to move between time-varying and time-constant partitions.
    bool _isEventTimeVariable;

public:

    // constructors, depending on whether you want trait rate or lambda/mu
    BranchEvent(double speciation, double lamshift, double extinction,
        double mushift, Node* x, Tree* tp, MbRandom* rp, double map,
        double scale);

    ~BranchEvent();

    void   setMapTime(double x);
    double getMapTime();

    void  setEventNode(Node* x);
    Node* getEventNode();

    void   setAbsoluteTime(double x);
    double getAbsoluteTime();

    void   setLamInit(double x);
    double getLamInit();

    void   setMuInit(double x);
    double getMuInit();

    void   setLamShift(double x);
    double getLamShift();

    void   setMuShift(double x);
    double getMuShift();

    void incrementMapPosition(double ink);
    void moveEventLocal();
    void moveEventGlobal();
    void setEventByMapPosition(double x);

    // Functions to set and manipulate OLD events:
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


inline void BranchEvent::setLamInit(double x)
{
    _lamInit = x;
}


inline double BranchEvent::getLamInit()
{
    return _lamInit;
}


inline void BranchEvent::setMuInit(double x)
{
    _muInit = x;
}


inline double BranchEvent::getMuInit()
{
    return _muInit;
}


inline void BranchEvent::setLamShift(double x)
{
    _lamShift = x;
}


inline double BranchEvent::getLamShift()
{
    return _lamShift;
}


inline void BranchEvent::setMuShift(double x)
{
    _muShift = x;
}


inline double BranchEvent::getMuShift()
{
    return _muShift;
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


inline bool BranchEvent::getIsEventTimeVariable()
{
    return _isEventTimeVariable;
}


inline void BranchEvent::setIsEventTimeVariable(bool x)
{
    _isEventTimeVariable = x;
}


#endif
