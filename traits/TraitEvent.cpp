/*
 *  TraitEvent.cpp
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
  *
 */

#include "TraitEvent.h"

 
#include <iostream>
#include <stdlib.h>

#include "TraitEvent.h"
#include "MbRandom.h"
#include "node.h"
 

// Made event time variable from the start with shift variable passed to constructor...


TraitBranchEvent::TraitBranchEvent(double beta, double shift, Node *x, Tree *tp, MbRandom* rp, double map, double scale){
	
	// Speciation and extinction from new event assumed constant in time...
	
	_betaInit = beta;
	_betaShift = shift;
 
	_absTime = tp->getAbsoluteTimeFromMapTime(map);
	
	nodeptr = x;
	mapTime = map;
	treePtr	= tp;
	ranPtr = rp;
	epsilon = scale;
	
	oldMapTime = map;
	oldNodePtr = x;
	
	_isEventTimeVariable = false; 
	
	
}



TraitBranchEvent::~TraitBranchEvent(void){
	
	
}

/* Oct 14 2012:
 Reversed sign of comparison operator-
 if maptime > event X maptime
	current event is LESS THAN or younger, than the focal event
 
 */

bool TraitBranchEvent::operator<(const TraitBranchEvent& a) const {
	
    if ( mapTime > a.mapTime )
        return true;
    return false;
}


/* this (along with incrementMapPosition) implements a uniform proposal mechanism
 with reflecting boundaries at all major tree boundaries: at the root, the proposal
 is bounced back onto the parent branch.
 At the tips OR at invalid nodes (defined in tree constructor)
 
 
 This mechanism of local move involves actually moving the position of the event
 If the proposed move is accepted, this is fine. If the proposed move is rejected,
 then the branchEvent will need to be reset to its original value. 
 IN practice, this should be done by:
 oldMapPosition = getMapTime();
 moveEventLocal();
 calculate likelihood etc
 if (accept)
 event has already been updated.
 if (reject)
 setEventByMapPosition(oldMapPosition)
 // should undo it.
 
 */



void TraitBranchEvent::moveEventLocal(void){
	
	oldNodePtr = nodeptr;
	oldMapTime = mapTime;
	
	double shift = ranPtr->uniformRv(0, epsilon) - 0.5*epsilon;
	//cout << "shifting event by " << shift << endl;
	
	incrementMapPosition(shift);
	
	// shouldn't actually need to update model attributes.
	// if node position shifts, it just brings its current attributes 
	// along, e.g., lambda etc.
	
	
}


void TraitBranchEvent::incrementMapPosition(double ink){
	
	double mapstart = nodeptr->getMapStart();
	double mapend = nodeptr->getMapEnd();
	
	
	double temp = getMapTime() + ink;
	
	// Map end is ROOTWARDS
	// Map start is TIPWARDS
	
	if (temp > mapend){
		// Goes to next most "ROOTWARDS" branch 
		// BUT if there is no next branch (because ancestor of curr branch is ROOT)
		//		go FORWARD down other descendant.
		
		ink = temp - mapend; // residual...
		
		if (getEventNode() == treePtr->getRoot()){
			cout << " Root problem in incrementMapPosition: should never get here" << endl;
			throw;
		}
		
		if (getEventNode()->getAnc() == treePtr->getRoot()){
			
			ink = -1*ink;		
			
			if (treePtr->getRoot()->getLfDesc() == getEventNode()){
				setEventNode(treePtr->getRoot()->getRtDesc());
				
			}else if (treePtr->getRoot()->getRtDesc() == getEventNode()){
				setEventNode(treePtr->getRoot()->getLfDesc());
				
			}else{
				cout << " error in incrementMapPosition()" << endl;
				throw;
			}
			setMapTime(getEventNode()->getMapEnd());	
			incrementMapPosition(ink);			
			
		}else{
			
			setEventNode(getEventNode()->getAnc());
			setMapTime(getEventNode()->getMapStart());		
			incrementMapPosition(ink);		
		}
		
		
		
		
	}else if (temp < mapstart){
		
		ink = temp - mapstart; // ink is now negative...
		
		if (getEventNode()->getLfDesc() == NULL && getEventNode()->getRtDesc() == NULL){
			// At terminal branch:
			// Reflect
			ink  = -1*ink; // make positive....
			setMapTime(mapstart);
			incrementMapPosition(ink);						
			
		}else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true && getEventNode()->getRtDesc()->getCanHoldEvent() == true){
			// both desc branches valid
			double ran = ranPtr->uniformRv();
			if (ran <= 0.5){
				// left branch
				setEventNode(getEventNode()->getLfDesc());
				setMapTime(getEventNode()->getMapEnd());
				
			}else{
				setEventNode(getEventNode()->getRtDesc());
				setMapTime(getEventNode()->getMapEnd());
			}
			incrementMapPosition(ink);			
		}else if(getEventNode()->getLfDesc()->getCanHoldEvent() == false && getEventNode()->getRtDesc()->getCanHoldEvent() == true){
			setEventNode(getEventNode()->getRtDesc());
			setMapTime(getEventNode()->getMapEnd());		
			incrementMapPosition(ink);	
			
		}else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true && getEventNode()->getRtDesc()->getCanHoldEvent()==false){
			
			setEventNode(getEventNode()->getLfDesc());
			setMapTime(getEventNode()->getMapEnd());		
			incrementMapPosition(ink);			
		}else if (getEventNode()->getLfDesc()->getCanHoldEvent() == false && getEventNode()->getRtDesc()->getCanHoldEvent() == false){
			// neither node can hold event: reflect
			ink  = -1*ink; // make positive....
			
			// Reflect:
			setMapTime(mapstart);
			incrementMapPosition(ink);			
			
		}else{
			cout << "Problem in incrementMapPosition()" << endl;
			throw;
		}	
		
		
	}else{
		
		setMapTime(temp);	
		
	}
	
	setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));	
	
	
}

/*
 
 This version has serious problems and does not move events correctly.
 The commented version below is the pre-October 13 2012 version.
 
 */

/*
void TraitBranchEvent::incrementMapPosition(double ink){
	
	double b_start = nodeptr->getMapStart();
	double b_end = nodeptr->getMapEnd();
	
	double temp = getMapTime() + ink;
	
	if (temp < b_start){
		if (getEventNode() == treePtr->getRoot()){
			cout << "Root problem: should never get here!" << endl;
		}
		
		ink = temp - b_start; 
		
		if (getEventNode()->getAnc() == treePtr->getRoot()){
			ink = -1*ink;
		}
		
		setEventNode(getEventNode()->getAnc());
		setMapTime(getEventNode()->getMapEnd());
		incrementMapPosition(ink);
		
	}else if (temp > b_end){
		
		ink = temp - b_end;
		
		if (getEventNode()->getLfDesc() == NULL && getEventNode()->getRtDesc() == NULL){
			// At terminal branch:
			// Reflect
			ink  = -1*ink;
			setMapTime(b_end);
			incrementMapPosition(ink);						
			
		}else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true && getEventNode()->getRtDesc()->getCanHoldEvent() == true){
			// both desc branches valid
			double ran = ranPtr->uniformRv();
			if (ran <= 0.5){
				// left branch
				setEventNode(getEventNode()->getLfDesc());
				setMapTime(getEventNode()->getMapStart());
				
			}else{
				setEventNode(getEventNode()->getRtDesc());
				setMapTime(getEventNode()->getMapStart());
			}
			incrementMapPosition(ink);			
		}else if(getEventNode()->getLfDesc()->getCanHoldEvent() == false && getEventNode()->getRtDesc()->getCanHoldEvent() == true){
			setEventNode(getEventNode()->getRtDesc());
			setMapTime(getEventNode()->getMapStart());		
			incrementMapPosition(ink);	
			
		}else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true && getEventNode()->getRtDesc()->getCanHoldEvent()==false){
			setEventNode(getEventNode()->getLfDesc());
			setMapTime(getEventNode()->getMapStart());		
			incrementMapPosition(ink);			
		}else if (getEventNode()->getLfDesc()->getCanHoldEvent() == false && getEventNode()->getRtDesc()->getCanHoldEvent() == false){
			// neither node can hold event: reflect
			ink  = -1*ink;
			setMapTime(b_end);
			incrementMapPosition(ink);			
			
		}else{
			cout << "Problem in incrementMapPosition()" << endl;
			throw;
		}
	}else{
		
		setMapTime(temp);
		
	}
	setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
	
}
*/


void TraitBranchEvent::moveEventGlobal(void){
	
	oldNodePtr = nodeptr;
	oldMapTime = mapTime;	
	
	double aa = treePtr->getRoot()->getMapStart();
	double bb = treePtr->getTotalMapLength();
	double position = ranPtr->uniformRv(aa, bb);
	setEventNode(treePtr->mapEventToTree(position));
	setMapTime(position);
	setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
	
}


// test
void TraitBranchEvent::setEventByMapPosition(double x){
	setEventNode(treePtr->mapEventToTree(x));
	setMapTime(x);
	setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));	
}


void TraitBranchEvent::revertOldMapPosition(void){
	
	nodeptr = oldNodePtr;
	mapTime = oldMapTime;
	setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
	
}
















