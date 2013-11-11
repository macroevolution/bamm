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
#include "Utilities.h"

class BranchEvent;
class comp_history;
 
typedef std::set<BranchEvent*, comp_history>::size_type EventSetSizeType;


class BranchHistory {

private:
	BranchEvent*					nodeEvent; // event describing focal node
	BranchEvent*					ancestralNodeEvent; // event describing ancestor...

	// set of all events on branch. Of length 0 if no events occurred on branch.
	// Also, if no events occur on branch, then entire branch is described by
	// the event referenced at nodeEvent;
	
	std::set<BranchEvent*,comp_history>	eventsOnBranch;  

	// New parameters, March 24 2012 for burst-speciation model
	double				_meanLambda;
	double				_meanMu;
			
	
	
		
public:
		
	BranchHistory(void);
	~BranchHistory(void);	
	
	
	BranchEvent*		getLastEvent(void); // get last event on branch
	
	// get last event from a reference event that occurred on branch:
	BranchEvent*		getLastEvent(BranchEvent* x); 
	BranchEvent*		getLastEvent(double ttime);
	BranchEvent*		getNextEvent(double ttime);
	int					getNumberOfEventsOnInterval(double t1, double t2);
	BranchEvent*		getEventByIndexPosition(int i);
	
	void				setNodeEvent(BranchEvent* x)				{ nodeEvent = x; }
	BranchEvent*		getNodeEvent(void)							{ return nodeEvent; }
	
	void				setAncestralNodeEvent(BranchEvent* x)		{ ancestralNodeEvent = x; }
	BranchEvent*		getAncestralNodeEvent(void)					{ return ancestralNodeEvent;  } 
	void				printBranchHistory(void);
	void				reversePrintBranchHistory(void);
	
	void				popEventOffBranchHistory(BranchEvent* x)	{ eventsOnBranch.erase(x); }
	void				addEventToBranchHistory(BranchEvent * x)	{ eventsOnBranch.insert(x); }
	int					getNumberOfBranchEvents(void)				{ return eventsOnBranch.size(); }
	
	
};	





#endif























