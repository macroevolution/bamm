 

#ifndef TRAITBRANCHHISTORY_H
#define TRAITBRANCHHISTORY_H

#include <set>

#include "TraitBranchEvent.h"
#include "Utilities.h"

class TraitBranchEvent;
class comp_history;

typedef std::set<TraitBranchEvent*, comp_history>::size_type EventSetSizeType;



class TraitBranchHistory {
	
private:
	TraitBranchEvent*					nodeEvent; // event describing focal node
	TraitBranchEvent*					ancestralNodeEvent; // event describing ancestor...
	
	// set of all events on branch. Of length 0 if no events occurred on branch.
	// Also, if no events occur on branch, then entire branch is described by
	// the event referenced at nodeEvent;
	
	std::set<TraitBranchEvent*,comp_history>	eventsOnBranch;  
	
	// New parameters, March 24 2012 for burst-speciation model
	double				_meanBeta;

public:
	
	TraitBranchHistory(void);
	~TraitBranchHistory(void);	
	
	
	TraitBranchEvent*		getLastEvent(void); // get last event on branch
	
	// get last event from a reference event that occurred on branch:
	TraitBranchEvent*		getLastEvent(TraitBranchEvent* x); 
	TraitBranchEvent*		getEventByIndexPosition(int i);
	
	void					setNodeEvent(TraitBranchEvent* x)			{ nodeEvent = x; }
	TraitBranchEvent*		getNodeEvent(void)							{ return nodeEvent; }
	
	void					setAncestralNodeEvent(TraitBranchEvent* x)	{ ancestralNodeEvent = x; }
	
	TraitBranchEvent*		getAncestralNodeEvent(void)					{ return ancestralNodeEvent;  } 
	void					printBranchHistory(void);
	void					reversePrintBranchHistory(void);
	
	void					popEventOffBranchHistory(TraitBranchEvent* x)	{ eventsOnBranch.erase(x); }
	void					addEventToBranchHistory(TraitBranchEvent * x)	{ eventsOnBranch.insert(x); }
	int						getNumberOfBranchEvents(void)					{ return eventsOnBranch.size(); }
	
	
};	

#endif























