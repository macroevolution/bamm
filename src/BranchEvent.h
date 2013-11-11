/*
 *  BranchEvent.h
 *  rateshift
 *
 *  Created by Dan Rabosky on 12/5/11.
	
 
 // Major update no more phenotypic evolution
 //	March 25 2012
 
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

class BranchEvent {
	
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
	double		_lamInit;		// Initial speciation rate at event
	double		_muInit;	    // Initial Mu rate at event
	double		_lamShift;		// magnitude & direction of speciation shift
	double		_muShift;		// magnitude & direction of mu shift
	
	/*****************/
	// New parameters June 12 2012
	//	to allow rjMCMC to move between time-varying and time-constant partitions.
	bool		_isEventTimeVariable;
	
	
	
public:
	

	// constructors, depending on whether you want trait rate or lambda/mu
	BranchEvent(double speciation, double lamshift, double extinction, double mushift, Node *x, Tree *tp, MbRandom* rp, double map, double scale);

	~BranchEvent(void);
	
	void		setMapTime(double x)		{ mapTime = x; }
	double		getMapTime(void)			{ return mapTime; }
	void		setEventNode(Node* x)		{ nodeptr = x; }
	Node*		getEventNode(void)			{ return nodeptr; }
 
	void		setAbsoluteTime(double x)		{ _absTime = x;		}
	double		getAbsoluteTime(void)			{ return _absTime;	}
	
	void		setLamInit(double x)		{ _lamInit = x;		}
	double		getLamInit(void)			{ return _lamInit;	}
	
	void		setMuInit(double x)				{ _muInit = x;			}
	double		getMuInit(void)					{ return _muInit;		}
	
	//void		setZparam(double x)			{ _zpar = x;		}
	//double	getZparam(void)				{ return _zpar;		}
	
	void		setLamShift(double x)		{ _lamShift = x;	}
	double		getLamShift(void)			{ return _lamShift;	}
	
	void		setMuShift(double x)		{ _muShift = x;		}
	double		getMuShift(void)			{ return _muShift;	}
	
	
	void		incrementMapPosition(double ink);
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
	bool		operator<(const BranchEvent& a) const;	
	
	// For time-varying rjMCMC:
	bool		getIsEventTimeVariable(void)	{ return _isEventTimeVariable;	}
	void		setIsEventTimeVariable(bool x)	{ _isEventTimeVariable = x;		}
	
	
};



#endif
