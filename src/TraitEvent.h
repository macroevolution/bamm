/*
 *  TraitEvent.h
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
  *
 */

/*
 *  event.h
 *  rateshift
 *
 *  Created by Dan Rabosky on 12/5/11.
 
 
 */

#ifndef TRAITEVENT_H
#define TRAITEVENT_H



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

class TraitBranchEvent {
	
private:
	
	double		mapTime;
	Node*		nodeptr;
	Tree*		treePtr;
	MbRandom*	ranPtr;
	
	// epsilon for local move
	double		epsilon;
	
	// Keep values of the old pointer and old maptime associated
	//	with the event for FAST reference 
	//	if a LOCAL proposal is rejected.
	
	Node*		oldNodePtr;
	double		oldMapTime;
	
	/******************/
	// New event parameters, March 25, 2012
	double		_absTime;	// real time, measured with t = 0 at root.
	double		_betaInit;	// initial beta value.
	double		_betaShift; // temporal shift parameter of trait evolution rate.
	/*****************/
	// New parameters June 12 2012
	//	to allow rjMCMC to move between time-varying and time-constant partitions.
	bool		_isEventTimeVariable;
	
	
	
public:
	
	
	// constructors, depending on whether you want trait rate or lambda/mu
	TraitBranchEvent(double beta, double shift, Node *x, Tree *tp, MbRandom* rp, double map, double scale);
	
	~TraitBranchEvent(void);
	
	void		setMapTime(double x)		{ mapTime = x; }
	double		getMapTime(void)			{ return mapTime; }
	void		setEventNode(Node* x)		{ nodeptr = x; }
	Node*		getEventNode(void)			{ return nodeptr; }
	
	void		setAbsoluteTime(double x)		{ _absTime = x;		}
	double		getAbsoluteTime(void)			{ return _absTime;	}
	
	void		setBetaInit(double x)		{	_betaInit = x;	}
	double		getBetaInit(void)			{	return _betaInit;	}
	
	void		setBetaShift(double x)		{  _betaShift = x;		}
	double		getBetaShift(void)			{	return _betaShift;	}
	
	void		incrementMapPosition(double ink);
	//void		incrementMapPosition(double ink, double starttime, Node * y);
	void		moveEventLocal(void);
	void		moveEventGlobal(void);
	void		setEventByMapPosition(double x);
	
	// functions to set and manipulate OLD events:
	void		setOldEventNode(Node* x)	{ oldNodePtr =  x; }
	Node*		getOldEventNode(void)		{ return oldNodePtr; }
	void		setOldMapTime(double x)		{ oldMapTime = x; }
	double		getOldMapTime(void)			{ return oldMapTime; }
	
	// revert to old map position using oldPtr and oldMapTime
	// this only works if you have changed the nodeptr and maptime
	//	relative to the values of oldNodePtr and oldMapTime;
	void		revertOldMapPosition(void);
	
	// overloading comparision operator:
	bool		operator<(const TraitBranchEvent& a) const;	
	
	// For time-varying rjMCMC:
	bool		getIsEventTimeVariable(void)	{ return _isEventTimeVariable;	}
	void		setIsEventTimeVariable(bool x)	{ _isEventTimeVariable = x;		}
	
	
};



#endif
