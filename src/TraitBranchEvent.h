/*
 *  TraitBranchEvent.h
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
  *
 */

#ifndef TRAITBRANCHEVENT_H
#define TRAITBRANCHEVENT_H

//Forward declarations:
class Tree;
class Node;
class MbRandom;

/*

 class TraitBranchEvent is the actual event.
 It contains:
 (1) the node associated with the event
 (2) the map position of the event
 (3) parameters associated with the event

 The initial event will always be the root node
 with a map position of 0.
 This event is immutable.

 */

class TraitBranchEvent
{

private:

    double mapTime;
    Node* nodeptr;
    Tree* treePtr;
    MbRandom* ranPtr;

    // epsilon for local move
    double epsilon;

    // Keep values of the old pointer and old maptime associated
    //  with the event for FAST reference
    //  if a LOCAL proposal is rejected.

    Node*  oldNodePtr;
    double oldMapTime;

    /******************/
    // New event parameters, March 25, 2012
    double _absTime; // real time, measured with t = 0 at root.
    double _betaInit;  // initial beta value.
    double _betaShift; // temporal shift parameter of trait evolution rate.

    /*****************/
    // New parameters June 12 2012
    // Allow rjMCMC to move between time-varying and time-constant partitions.
    bool _isEventTimeVariable;

public:

    // constructors, depending on whether you want trait rate or lambda/mu
    TraitBranchEvent(double beta, double shift, Node* x, Tree* tp, MbRandom* rp,
                     double map, double scale);

    ~TraitBranchEvent();

    void   setMapTime(double x);
    double getMapTime();

    void  setEventNode(Node* x);
    Node* getEventNode();

    void   setAbsoluteTime(double x);
    double getAbsoluteTime();

    void   setBetaInit(double x);
    double getBetaInit();

    void   setBetaShift(double x);
    double getBetaShift();

    void incrementMapPosition(double ink);
    void moveEventLocal();
    void moveEventGlobal();
    void setEventByMapPosition(double x);

    // functions to set and manipulate OLD events:
    void  setOldEventNode(Node* x);
    Node* getOldEventNode();

    void   setOldMapTime(double x);
    double getOldMapTime();

    // revert to old map position using oldPtr and oldMapTime
    // this only works if you have changed the nodeptr and maptime
    // relative to the values of oldNodePtr and oldMapTime;
    void revertOldMapPosition();

    // overloading comparision operator:
    bool operator<(const TraitBranchEvent& a) const;

    // For time-varying rjMCMC:
    bool getIsEventTimeVariable();
    void setIsEventTimeVariable(bool x);
};


inline void TraitBranchEvent::setMapTime(double x)
{
    mapTime = x;
}


inline double TraitBranchEvent::getMapTime()
{
    return mapTime;
}


inline void TraitBranchEvent::setEventNode(Node* x)
{
    nodeptr = x;
}


inline Node* TraitBranchEvent::getEventNode()
{
    return nodeptr;
}


inline void TraitBranchEvent::setAbsoluteTime(double x)
{
    _absTime = x;
}


inline double TraitBranchEvent::getAbsoluteTime()
{
    return _absTime;
}


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


inline void TraitBranchEvent::setOldEventNode(Node* x)
{
    oldNodePtr =  x;
}


inline Node* TraitBranchEvent::getOldEventNode()
{
    return oldNodePtr;
}


inline void TraitBranchEvent::setOldMapTime(double x)
{
    oldMapTime = x;
}


inline double TraitBranchEvent::getOldMapTime()
{
    return oldMapTime;
}


inline bool TraitBranchEvent::getIsEventTimeVariable()
{
    return _isEventTimeVariable;
}


inline void TraitBranchEvent::setIsEventTimeVariable(bool x)
{
    _isEventTimeVariable = x;
}


#endif
