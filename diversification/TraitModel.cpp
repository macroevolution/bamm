/*
 *  TraitModel.cpp
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
 *
 */

#include "TraitModel.h"


#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <fstream>
#include <limits>

#include "MbRandom.h"
#include "node.h"
#include "TraitBranchHistory.h"
#include "Settings.h"
#include "Utilities.h"

//#define DEBUG_LH

double TraitModel::mhColdness = 1.0;

TraitModel::TraitModel(MbRandom * ranptr, Tree * tp, Settings * sp){
	
	_lastLH = 0.0;
	
	cout << endl << "Initializing model object...." << endl;
	
	// reduce weird autocorrelation of values at start by calling RNG a few times...
	for (int i =0; i<100; i++)
		double x = ranptr->uniformRv();
	
	ran = ranptr;
	treePtr = tp;
	sttings = sp;
	
	//cout << "model ML: " << treePtr->getTotalMapLength() << endl;
	
	double bl = treePtr->getTotalMapLength(); // total map length (required to set priors)
	
	// Set parameter values for model object, including priors etc.
	
	gen = 0;
	
	/*	******************	*/ 
	/* MCMC proposal/tuning parameters  */
	/*
	 I am keeping all of these parameters with apparent double-definitions:
	 They are stored both in Settings (sttings)
	 but also as private data for class Model
	 
	 The rationale for this is that I may keep the class Settings as a 
	 one-way class that holds INITIAL settings
	 But may incorporate some auto-tuning of parameters to make the MCMC 
	 more efficient, and thus the within-Model parameters will be updated
	 but not the Settings parameters. 
	 
	 This could all change, though.
	 */

	_updateBetaScale = sttings->getUpdateBetaScale();
	_updateBetaShiftScale = sttings->getUpdateBetaShiftScale();
	_updateNodeStateScale = sttings->getUpdateNodeStateScale();
	
	// initial values

	
	/*
	 Scale parameter for step size is determined dynamically by average waiting time 
	 between successive speciation events on tree.
	 A value of _MeanSpeciationLengthFraction = 0.1 means that step sizes will be chosen from
	 a uniform (-x, x) distribution centered on the current location, where 
	 x is 0.10 times the mean waiting time (or total tree length / # speciation events
	 */
	
	_scale = sttings->getMeanSpeciationLengthFraction() * (treePtr->getTreeLength()); // scale for event moves on tree.
	_scale /= (double)(treePtr->getNumberTips());
	
	//cout << "SCALE: " << _scale << endl;
	
	_updateEventRateScale = sttings->getUpdateEventRateScale();
	_localGlobalMoveRatio = sttings->getLocalGlobalMoveRatio();	// For Poisson process
	_targetNumber = sttings->getTargetNumberOfEvents();
	
	// Keep poisson rate prior such that expected event rate will generate 1 event on the tree.		
	//poissonRatePrior = bl; // important to keep # of events low.
	// setting this equal to bl leads to expected number of 1 event
	poissonRatePrior = bl / _targetNumber; // should lead to TARGET number of events from prior...
	
	// Keep priors in settings...
	setMinMaxTraitPriors();
	
	eventLambda = (1/bl); // event rate, initialized to generate expected number of 1 event
 
	//cout << bl << "\t" << eventLambda << endl;
	
	/*	*********************	*/ 
	/* Other parameters & tracking variables*/ 
	
	acceptCount = 0;
	rejectCount = 0;
	acceptLast = -1;
	
	//set up event at root node:
	double startTime = 0;
	
	TraitBranchEvent* x =  new TraitBranchEvent((double)sttings->getBetaInit(), sttings->getBetaShiftInit(), treePtr->getRoot(), treePtr, ran, startTime, _scale);
	rootEvent = x;
	lastEventModified = x;		
	
	cout << "Root beta: " << x->getBetaInit() << "\t" << sttings->getBetaInit() << "\tShift: ";
	cout << x->getBetaShift() << endl;
	
	// set NodeEvent of root node equal to the rootEvent:
	tp->getRoot()->getTraitBranchHistory()->setNodeEvent(rootEvent);	
	
	//initializing all branch histories to equal the root event:
	forwardSetBranchHistories(rootEvent);
	
	// 
	treePtr->setMeanBranchTraitRates();
	
	setCurrLnLTraits(computeLikelihoodTraits());
	
 
	cout << "Model object successfully initialized." << endl;
	cout << "Initial log-likelihood: " << getCurrLnLTraits() << endl << endl;
	if (sttings->getSampleFromPriorOnly()){
		cout << "\tNote that you have chosen to sample from prior only." << endl;
	}
	
}


TraitModel::~TraitModel(void){
	
	for (std::set<TraitBranchEvent*>::iterator it = eventCollection.begin(); it != eventCollection.end(); it++)
		delete (*it);
}

/*
 Adds event to tree based on reference map value
 -adds to branch history set
 -inserts into Model::eventCollection
 
 
 
 */


void TraitModel::addEventToTree(double x){
	
	
	// For now, the rates of speciation and extinction are set to whatever they should be based
	// on the ancestralNodeEvent
	Node * xnode = treePtr->mapEventToTree(x);
	double atime = treePtr->getAbsoluteTimeFromMapTime(x);
	TraitBranchHistory * bh = xnode->getTraitBranchHistory();
	TraitBranchEvent * be = bh->getAncestralNodeEvent();
	
	double elapsed = atime - be->getAbsoluteTime();
	double newbeta = be->getBetaInit() * exp( elapsed * be->getBetaShift());
 
	// End calculations:: now create event
	
	TraitBranchEvent* newEvent = new TraitBranchEvent(newbeta,  be->getBetaShift(), treePtr->mapEventToTree(x), treePtr, ran, x, _scale);
	
	// add the event to the branch history. 
	//	ALWAYS done after event is added to tree.
	newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(newEvent);
	
	eventCollection.insert(newEvent);	
	
	// Event is now inserted into branch history: 
	//	however, branch histories must be updated.
	
	forwardSetBranchHistories(newEvent);
	
	treePtr->setMeanBranchTraitRates();
	
	// Addition June17 2012
	lastEventModified = newEvent;	
	
}


/*
 Adds event to tree based on uniform RV
 -adds to branch history set
 -inserts into Model::eventCollection
 
 
 
 */


void TraitModel::addEventToTree(void){
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\tenter add event ot tree\t" << computeLikelihoodTraits() << endl;
	
#endif			
	
	double aa = treePtr->getRoot()->getMapStart();
	double bb = treePtr->getTotalMapLength();
	double x = ran->uniformRv(aa, bb);
	
	
	// For now, the rates of speciation and extinction are set to whatever they should be based
	// on the ancestralNodeEvent
	Node * xnode = treePtr->mapEventToTree(x);
	double atime = treePtr->getAbsoluteTimeFromMapTime(x);
	TraitBranchHistory * bh = xnode->getTraitBranchHistory();
	TraitBranchEvent * be = bh->getAncestralNodeEvent();
	
	double elapsed = atime - be->getAbsoluteTime();
 
	double newbeta = be->getBetaInit() * exp( elapsed * be->getBetaShift());
	
	// End calculations:: now create event
	
	TraitBranchEvent* newEvent = new TraitBranchEvent(newbeta, be->getBetaShift(), treePtr->mapEventToTree(x), treePtr, ran, x, _scale);
	
	newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(newEvent);
	
	eventCollection.insert(newEvent);
	
	// Event is now inserted into branch history: 
	//	however, branch histories must be updated.
	
	forwardSetBranchHistories(newEvent);
	treePtr->setMeanBranchTraitRates();
	
	// Addition June17 2012
	lastEventModified = newEvent;
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\t Leave addEventToTree\t" << computeLikelihoodTraits() << endl;
	
#endif	
	
}


// This function for adding event with beta...

void TraitModel::addEventToTreeWithSetBeta(double beta, double bshift){
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\tenter add event ot tree\t" << computeLikelihoodTraits() << endl;
	
#endif			
	
	double aa = treePtr->getRoot()->getMapStart();
	double bb = treePtr->getTotalMapLength();
	double x = ran->uniformRv(aa, bb);
	
	
	// For now, the rates of speciation and extinction are set to whatever they should be based
	// on the ancestralNodeEvent
	Node * xnode = treePtr->mapEventToTree(x);
	double atime = treePtr->getAbsoluteTimeFromMapTime(x);
	TraitBranchHistory * bh = xnode->getTraitBranchHistory();
	//TraitBranchEvent * be = bh->getAncestralNodeEvent();
	
	//double elapsed = atime - be->getAbsoluteTime();
 
	// End calculations:: now create event
	
	TraitBranchEvent* newEvent = new TraitBranchEvent(beta, bshift, treePtr->mapEventToTree(x), treePtr, ran, x, _scale);
	
	newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(newEvent);
	
	eventCollection.insert(newEvent);
	
	// Event is now inserted into branch history: 
	//	however, branch histories must be updated.
	
	forwardSetBranchHistories(newEvent);
	treePtr->setMeanBranchTraitRates();
	
	// Addition June17 2012
	lastEventModified = newEvent;
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\t Leave addEventToTree\t" << computeLikelihoodTraits() << endl;
	
#endif		
}



void TraitModel::printEvents(void){
	
	// for each event:
	//	print:	maptime
	//			nodeptr
	//			
	int n_events = eventCollection.size();
	cout << "N_events: " << n_events << endl;
	int counter = 1;
	for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin(); i != eventCollection.end(); i++){
		cout << "event " << counter++ << "\tAddress: " << (*i) << "\t";
		cout << (*i)->getMapTime() << "\tNode: " << (*i)->getEventNode() << endl << endl;
	}
	
	
}

TraitBranchEvent* TraitModel::chooseEventAtRandom(void){
	
	int n_events = eventCollection.size();
	if (n_events == 0){
		return NULL;
		//should ultimately throw exception here.
		
	} else{
		int ctr = 0;
		double xx = ran->uniformRv();
		int chosen = (int)(xx * (double)n_events);		
		
		std::set<TraitBranchEvent*>::iterator sit = eventCollection.begin();
		
		for (int i = 0; i < chosen; i++){
			sit++;
			ctr++;
		}
		return (*sit);
	}
	
	
}





/*
 void Model::eventLocalMove(void)
 
 IF events are on tree:
 choose event at random
 move locally and forward set branch histories etc.
 should also store previous event information to revert to previous
 
 */

void TraitModel::eventLocalMove(void){

#ifdef DEBUG_LH
	cout << getGeneration() << "\tenter eventLocalMove\t" << computeLikelihoodTraits() << endl;
	
#endif
	
	
	
	if (getNumberOfEvents() > 0){
		
		// the event to be moved
		TraitBranchEvent* chosenEvent = chooseEventAtRandom();
		
		// corresponding node defining branch on which event occurs
		Node* theEventNode = chosenEvent->getEventNode();
		
		// this is the event preceding the chosen event: histories should be set forward from here..
		TraitBranchEvent* previousEvent = chosenEvent->getEventNode()->getTraitBranchHistory()->getLastEvent(chosenEvent);
		
		// set this history variable in case move is rejected
		lastEventModified = chosenEvent;
		
		chosenEvent->getEventNode()->getTraitBranchHistory()->popEventOffBranchHistory(chosenEvent);
		chosenEvent->moveEventLocal(); // move event
		chosenEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(chosenEvent);
		
		// Get last event from the theEventNode
		//		forward set its history
		//	Then go to the "moved" event and forward set its history
		
		
		//BranchEvent* newLastEvent = getLastEvent(theEventNode);
		
		
		forwardSetBranchHistories(previousEvent);
		forwardSetBranchHistories(chosenEvent);
		
	}
	// else no events to move
	//cout << "leave localMove" << endl;
	
	treePtr->setMeanBranchTraitRates();

#ifdef DEBUG_LH
	cout << getGeneration() << "\tleave eventLocalMove\t" << computeLikelihoodTraits() << endl;
	
#endif	
	
}

void TraitModel::eventGlobalMove(void){

#ifdef DEBUG_LH
	cout << getGeneration() << "\tenter eventGlobalMove\t" << computeLikelihoodTraits() << endl;
	
#endif		
	
	if (getNumberOfEvents() > 0){
		TraitBranchEvent* chosenEvent = chooseEventAtRandom();
		
		// this is the event preceding the chosen event: histories should be set forward from here..
		TraitBranchEvent* previousEvent = chosenEvent->getEventNode()->getTraitBranchHistory()->getLastEvent(chosenEvent);
		
		//cout << "EGM: moving " << chosenEvent << "\tLastEvent: " << previousEvent << endl;
		
		Node* theEventNode = chosenEvent->getEventNode();
		
		// private variable
		lastEventModified = chosenEvent;
		
		chosenEvent->getEventNode()->getTraitBranchHistory()->popEventOffBranchHistory(chosenEvent);
		chosenEvent->moveEventGlobal(); // move event
		chosenEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(chosenEvent);
		
		// Get last event from the theEventNode
		//		forward set its history
		//	Then go to the "moved" event and forward set its history
		//BranchEvent* newLastEvent = getLastEvent(theEventNode);
		//cout << "EGM2: moving " << chosenEvent << "\tLastEvent: " << previousEvent <<  "\tBLE: " << newLastEvent << endl;
		
		
		forwardSetBranchHistories(previousEvent);
		forwardSetBranchHistories(chosenEvent);
		
	}
	//cout << "leave globalMove" << endl;
	
	treePtr->setMeanBranchTraitRates();

#ifdef DEBUG_LH
	cout << getGeneration() << "\tleave eventGlobalMove\t" << computeLikelihoodTraits() << endl;
	
#endif		
	
}

// used to reset position of event if move is rejected

void TraitModel::revertMovedEventToPrevious(void){

#ifdef DEBUG_LH
	cout << getGeneration() << "\tenter revertMovedEvntToPrv\t" << computeLikelihoodTraits() << endl;
	
#endif		
	
	double startLH = getCurrLnLTraits();
	
	// Get LAST EVENT from position of event to be removed:
	
	TraitBranchEvent * newLastEvent = lastEventModified->getEventNode()->getTraitBranchHistory()->getLastEvent(lastEventModified);
	
	//BranchEvent * newLastEvent = getLastEvent(lastEventModified);	
	
	// pop event off its new position
	lastEventModified->getEventNode()->getTraitBranchHistory()->popEventOffBranchHistory(lastEventModified);
	
	// Reset nodeptr:
	// Reset mapTime:
	lastEventModified->revertOldMapPosition();
	
	// Now: reset forward from lastEventModified (new position)
	//	and from newLastEvent, which holds 'last' event before old position
	
	lastEventModified->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(lastEventModified);
	
	// Forward set from new position
	forwardSetBranchHistories(newLastEvent);
	
	// forward set from event immediately rootwards from previous position:
	forwardSetBranchHistories(lastEventModified);
	
	
	
	// Set lastEventModified to NULL,
	//	because it has already been reset.
	//	Future implementations should check whether this is NULL
	//	before attempting to use it to set event
	
	lastEventModified = NULL;
	
	// Reset speciaton-extinction on branches
	treePtr->setMeanBranchTraitRates();
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\tleave revertMovedEvntToPrv\t" << computeLikelihoodTraits() << endl;
	
#endif	
	
}



// Recursively count the number of events in the branch histories
int TraitModel::countEventsInBranchHistory(Node* p){
	int count = 0;
	count += p->getTraitBranchHistory()->getNumberOfBranchEvents();
	if (p->getLfDesc() != NULL){
		count += countEventsInBranchHistory(p->getLfDesc());
	}
	if (p->getRtDesc() != NULL){
		count += countEventsInBranchHistory(p->getRtDesc());
	}
	
	return count;
}

/*
 
 Deletes an event from tree.
 
 */

void TraitModel::deleteEventFromTree(TraitBranchEvent * be){
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\tenter deleteEventFromTree(be)\t" << computeLikelihoodTraits() << endl;
	
#endif			
	
	if (be == rootEvent){
		cout << "Can't delete root event" << endl;
		exit(1);
	}else{
		// erase from branch history:
		Node* currNode = (be)->getEventNode();
		
		//get event downstream of i
		TraitBranchEvent * newLastEvent = currNode->getTraitBranchHistory()->getLastEvent(be);
		
		lastDeletedEventMapTime = (be)->getMapTime();
		
		_lastDeletedEventBetaInit = (be)->getBetaInit();
		_lastDeletedEventBetaShift = (be)->getBetaShift();
 
		
		currNode->getTraitBranchHistory()->popEventOffBranchHistory((be));
		
		// delete from global node set
		delete (be);				
		//cout << "deleted..." << endl;
		
		eventCollection.erase(be);
		
		forwardSetBranchHistories(newLastEvent);
		
	}
	
	
	treePtr->setMeanBranchTraitRates();
#ifdef DEBUG_LH
	cout << getGeneration() << "\tleave deleteEventFromTree(be)\t" << computeLikelihoodTraits() << endl;
	
#endif		

}




void TraitModel::deleteRandomEventFromTree(void){
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\temter deleteRandomEventFromTree(void)\t" << computeLikelihoodTraits() << endl;
	
#endif		
	
	
	//cout << endl << endl << "START Delete: " << endl;
	//printBranchHistories(treePtr->getRoot());
	
	// can only delete event if more than root node present.
	int n_events = eventCollection.size();
	
	if (eventCollection.size() > 0){
		int counter = 0;
		double xx = ran->uniformRv();
		int chosen = (int)(xx * (double)n_events);	
		
		for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin(); i != eventCollection.end(); i++){
			if (counter++ == chosen){
				
				// erase from branch history:
				Node* currNode = (*i)->getEventNode();
				
				//get event downstream of i
				TraitBranchEvent * newLastEvent = currNode->getTraitBranchHistory()->getLastEvent((*i));
				
				lastDeletedEventMapTime = (*i)->getMapTime();
				//lastDeletedEventBeta = (*i)->getBeta();
				
				_lastDeletedEventBetaInit = (*i)->getBetaInit();
				_lastDeletedEventBetaShift = (*i)->getBetaShift();
 
				
				currNode->getTraitBranchHistory()->popEventOffBranchHistory((*i));
				
				
				// delete from global node set
				delete (*i);				
				//cout << "deleted..." << endl;
				
				eventCollection.erase(i);
				
				//cout << "erased ... " << endl;
				
				
				// reset forward history from last event:
				//BranchEvent* lastEvent = getLastEvent(currNode);
				//forwardSetBranchHistories(lastEvent);
				forwardSetBranchHistories(newLastEvent);
				
				//cout << "forward set..." << endl;
				
				// is this correctly setting branch histories??
				
			}
		}
	}
	treePtr->setMeanBranchTraitRates();
#ifdef DEBUG_LH
	cout << getGeneration() << "\tLeave deleteRandomEventFromTree(void)\t" << computeLikelihoodTraits() << endl;
	
#endif		
}


// Valid, March 23

void TraitModel::restoreLastDeletedEvent(void){

#ifdef DEBUG_LH
	cout << getGeneration() << "\tEnter restoreLastDeletedEvent\t" << computeLikelihoodTraits() << endl;
	
#endif		
	
	// Constructor for traitEvolution model:
	
	//BranchEvent* newEvent = new BranchEvent(lastDeletedEventBeta, treePtr->mapEventToTree(lastDeletedEventMapTime), treePtr, ran, lastDeletedEventMapTime);
	
	// Use constructor for speciation and extinction
	
	TraitBranchEvent * newEvent = new TraitBranchEvent((double)0.0, (double)0.0, treePtr->mapEventToTree(lastDeletedEventMapTime), treePtr, ran, lastDeletedEventMapTime, _scale);
	
	newEvent->setBetaInit(_lastDeletedEventBetaInit);
	newEvent->setBetaShift(_lastDeletedEventBetaShift);
	
	// add the event to the branch history. 
	//	ALWAYS done after event is added to tree.
	newEvent->getEventNode()->getTraitBranchHistory()->addEventToBranchHistory(newEvent);
	
	eventCollection.insert(newEvent);	
	
	// Event is now inserted into branch history: 
	//	however, branch histories must be updated.
	
	forwardSetBranchHistories(newEvent);
	
	treePtr->setMeanBranchTraitRates();
	//setCurrLnLTraits(computeLikelihoodTraits());
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\tLeave restoreLastDeletedEvent\t" << computeLikelihoodTraits() << endl;
	
#endif		
}

/*
 
 
 
 
 
 june20 through here update trait stuff.
 must go through and check for (remove)
 all references to getBranchHistory and related suff
 
 
 
 
 
 */ 


void TraitModel::changeNumberOfEventsMH(void){
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\tEnter chnageNumberEvetnsMH\t" << computeLikelihoodTraits() << endl;
	
#endif		
	
	// Propose gains & losses equally if not on boundary (n = 0) events:
	
	bool gain = (ran->uniformRv() <= 0.5);		
	if (eventCollection.size() == 0){
		// set event to gain IF on boundary
		gain = true;
	} 
	
	// now to adjust acceptance ratios:
	
	if (gain){
		
		
		double qratio = 1.0;
		if (eventCollection.size() == 0){
			// no events on tree
			// can only propose gains.
			qratio = 0.5;
		}else{
			// DO nothing.
		}
		
#ifdef NO_DATA
		double likTraits = 0.0;
		double PropLnLik = likTraits;
		
#else
		treePtr->setMeanBranchTraitRates();
		double likTraits = computeLikelihoodTraits();
		double PropLnLik = likTraits;	
		
#endif			
		double logHR = log(treePtr->getTotalMapLength() * eventLambda * qratio);
		logHR -= log( (double)(eventCollection.size()+1));						
		
		double likeRatio = PropLnLik - getCurrLnLTraits();
		
		logHR += likeRatio;		
	
#ifdef DEBUG_LH
		cout << getGeneration() << "\tLogHR: " << logHR << "\tInitL: " << getCurrLnLTraits();
		cout << "\tNewL: "<< PropLnLik << endl;
#endif		
		
		
		bool acceptMove = false;
		if (std::isinf(likeRatio) ){
			
		}else {
			acceptMove = acceptMetropolisHastings(logHR);
		} 
		
 
		
		if (acceptMove){
			//cout << "gaining event in changeNumberOfEventsMH " << eventCollection.size() << "\t" << getCurrLnLTraits() << endl;
			//cout << "gainaccept" << computeLikelihoodBranches()  << endl;					
			
			addEventToTree();
			//cout << "Calliing isValid from ChangeNumberEvents::gain" << endl;
			
			bool isValidConfig = isEventConfigurationValid(lastEventModified);
			
			if (isValidConfig){
#ifdef DEBUG_LH
				cout << getGeneration() << "AcceptingGain: " << endl;
#endif				
				setCurrLnLTraits(likTraits);
				// Update accept/reject statistics
				acceptCount++;
				acceptLast = 1;
 			}else{
				
#ifdef DEBUG_LH
				cout << getGeneration() << "Rejecting Gain : invalid config: " << endl;
#endif		
				
				// Need to get rid of event that was just gained...
				//cout << "Invalid event config from addEventToTree - deleting." << endl;
				deleteEventFromTree(lastEventModified);
				treePtr->setMeanBranchTraitRates();
				rejectCount++;
				acceptLast = 0;			
 			}
			
			
		}else{
			
#ifdef DEBUG_LH
			cout << getGeneration() << "Rejecting Gain : MH reject: " << endl;
#endif		
			
			//cout << "gainreject" << computeLikelihoodBranches() << endl;
			rejectCount++;
			acceptLast = 0;
		}
		//setCurrLnLTraits(computeLikelihoodTraits());
		
		
		
		
	}else{
		//cout << "loss: initial LH: " << computeLikelihoodBranches() << endl;
		deleteRandomEventFromTree(); 
		//cout << "loss: LH after deleteRandomEvent" << computeLikelihoodBranches() << endl;
		
#ifdef NO_DATA
		double likTraits = 0.0;
		double PropLnLik = likTraits;
		
#else
		
		treePtr->setMeanBranchTraitRates();
		
		double likTraits = computeLikelihoodTraits();
		double PropLnLik = likTraits;	
		
#endif				
		
		double qratio = 1.0; // if loss, can only be qratio of 1.0
		
		if (eventCollection.size()  == 1){
			qratio = 2.0;
		} 
		
		// Note here that ratio is inverted from above...
		
		double logHR = log( (double)(eventCollection.size()+1) * qratio);	
		
		logHR -= log(treePtr->getTotalMapLength()* eventLambda );
		
		
		
		
		// Hasting ratio is just inverse for loss event
		// but likelihood ratio must be taken into account:
		
		double likeRatio = PropLnLik - getCurrLnLTraits();
		
		logHR += likeRatio;
		
		bool acceptMove = false;
		
		if (std::isinf(likeRatio) ){
			
		}else{
			acceptMove = acceptMetropolisHastings(logHR);
		} 
		
 
		
		//cout << "loss: " << acceptMove << "\t" << PropLnLik << "\tLT " << getCurrLnLTraits() + getCurrLnLBranches() << endl;
		
		if (acceptMove){
#ifdef DEBUG_LH
			cout << getGeneration() << "AcceptDelete: " << endl;
#endif
			
			//cout << "loss accept, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << endl;
			setCurrLnLTraits(likTraits);
			
			acceptCount++;
			acceptLast = 1;
			
			
		}else{
			
#ifdef DEBUG_LH
			cout << getGeneration() << "Rejecting Delete: " << endl;
#endif			
			restoreLastDeletedEvent();
			
			// Trait evolution rates on branches automatically updated after restoreLastDeletedEvent()
			setCurrLnLTraits(computeLikelihoodTraits());
			
			// LIkelihood should still need updated:
			
#ifdef DEBUG_LH
			
			double L1 = getCurrLnLTraits();
			double L2 = computeLikelihoodTraits();
			double delta = L1 - L2;
			
			if ((delta < -0.01) || (delta > 0.01)){
				cout << "ChangeNumberEventsMH:";
				cout << "problem at " << getGeneration() << "\tL1: " << getCurrLnLTraits();
				cout << "\tL2: " << computeLikelihoodTraits() << endl;
				exit(1);
			} 
#endif			
			
			
			//cout << "loss reject restored, LH: " << computeLikelihoodBranches() << "\tlikBranches" << likBranches << endl;		
			rejectCount++;
			acceptLast = 0;
			
			
		}	
		//setCurrLnLTraits(computeLikelihoodTraits());
	}
	
	// For now, assume Likelihood ratio is always 1 here
	// event gain or loss just drops new event on tree
	// but rates before and after event are identical until updated.
	
#ifdef DEBUG_LH
	cout << getGeneration() << "\tLeaving chnageNumberEvetnsMH\t" << computeLikelihoodTraits() << endl;
	
#endif			
	
	incrementGeneration();
	//cout << "leave changeNumberEvents" << endl;

}

void TraitModel::moveEventMH(void){
 
	
	if (eventCollection.size() > 0){
		
		double localMoveProb = _localGlobalMoveRatio / (1 + _localGlobalMoveRatio);
		
		bool isLocalMove = (ran->uniformRv() <= localMoveProb);
		//cout << "is local: " << isLocalMove << endl;
		
		if (isLocalMove){
			// Local move, with event drawn at random
			eventLocalMove();
			
#ifdef NO_DATA
			double likTraits = 0;
			double PropLnLik = likTraits;
#else
			
			treePtr->setMeanBranchTraitRates();
			
			double likTraits = computeLikelihoodTraits();
			double PropLnLik = likTraits;
			
			
			
#endif		
			
			double likeRatio = PropLnLik - getCurrLnLTraits();
			double logHR = likeRatio;
			
			// No longer protecting this as const bool
			
			bool acceptMove = false;
			bool isValid = false;
			//cout << "calling isValid from moveEventMH::local" << endl;
			isValid = isEventConfigurationValid(lastEventModified);
			
			if (std::isinf(likeRatio) ){
				
			}else if (isValid){
				acceptMove = acceptMetropolisHastings(logHR);
			}else{
				//cout << "Invalid event configuration from LocalMove" << endl;
			}		
			
			//const bool acceptMove = acceptMetropolisHastings(logHR);			
			
			if (acceptMove == true){
				setCurrLnLTraits(likTraits);
	
				acceptCount++;
				acceptLast = 1;
				
			}else{
				// revert to previous state
				revertMovedEventToPrevious();
				
				treePtr->setMeanBranchTraitRates();
  
				rejectCount++;
				acceptLast = 0;
			
				
#ifdef DEBUG_LH
			
				double L1 = getCurrLnLTraits();
				double L2 = computeLikelihoodTraits();
				double delta = L1 - L2;
				
				if ((delta < -0.01) || (delta > 0.01)){
					cout << "MoveEventLocal:";
					cout << "problem at " << getGeneration() << "\tL1: " << getCurrLnLTraits();
					cout << "\tL2: " << computeLikelihoodTraits() << "\tStart: " << startLH << endl;
					exit(1);
				} 
	 
#endif			
				
			
			}	
			
			
			
		}else{
			// global move, event drawn at random
			eventGlobalMove();
			//cout << "successful global move" << endl;
#ifdef NO_DATA
			double likTraits = 0;
			double PropLnLik = likTraits;
#else
			
			treePtr->setMeanBranchTraitRates();
			
			double likTraits = computeLikelihoodTraits();
			double PropLnLik = likTraits;
			
#endif		
			
			double likeRatio = PropLnLik - getCurrLnLTraits();
			double logHR = likeRatio;
	
			 //const bool acceptMove = acceptMetropolisHastings(logHR);
		 
			bool acceptMove = false;
			bool isValid = false;
			//cout << "calling isValid from moveEventMH::global" << endl;
			isValid = isEventConfigurationValid(lastEventModified);
			
			if (std::isinf(likeRatio) ){
				
			}else if (isValid){
				acceptMove = acceptMetropolisHastings(logHR);
			}else{
				//cout << "Invalid event configuration from GlobalMove" << endl;
			}				
			
			if (acceptMove == true){
 				setCurrLnLTraits(likTraits);
				acceptCount++;
				acceptLast = 1;
				
			}else{
				// revert to previous state
				revertMovedEventToPrevious();	
				
				treePtr->setMeanBranchTraitRates();
				rejectCount++;
				acceptLast = 0;
				
				
#ifdef DEBUG_LH
				
				double L1 = getCurrLnLTraits();
				double L2 = computeLikelihoodTraits();
				double delta = L1 - L2;
				
				if ((delta < -0.01) || (delta > 0.01)){
					cout << "MoveEventGlobal:";
					cout << "problem at " << getGeneration() << "\tL1: " << getCurrLnLTraits();
					cout << "\tL2: " << computeLikelihoodTraits() << endl;
					exit(1);
				} 
#endif			
				
				
				
			}	
			
			
		}
		
	}else{
		// consider proposal rejected (can't move nonexistent event)
		rejectCount++;
		acceptLast = 0;
	}
	
	
	incrementGeneration();
	
}




/* June 12 2012
 Select an event at random.
 If partition is time-constant
 flip state to time-variable
 If partition is time-variable
 flip state to time-constant
 
 */
void TraitModel::updateTimeVariablePartitionsMH(void){
	
	int n_events = eventCollection.size() + 1;
	int toUpdate = ran->sampleInteger(0, eventCollection.size());
	TraitBranchEvent* be = rootEvent;
	
	if (toUpdate > 0){
		std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
		for (int i = 1; i < toUpdate; i++){
			myIt++;		
		}
		
		be = (*myIt);
	}else{
		// event remains as root event-
	}	
	
	if (be->getIsEventTimeVariable()){
		
		
		
		
	}else if (!be->getIsEventTimeVariable()){
		
		
		
	}else{
		// Should not be able to get here:
		cout << "Invalid _isEventTimeVariable in Model::UpdateTimeVariablePartitionsMH" << endl;
		throw;
	}
	
	
}



/*
 
 Metropolis-Hastings step to update Poisson event rate.
 Note that changing this rate does not affect the likelihood,
 so the priors and qratio determine acceptance rate.
 
 */
void TraitModel::updateEventRateMH(void){
	
	//cout << "Entering update event rate" << endl;
	
	double oldEventRate = getEventRate();
	double cterm = exp( _updateEventRateScale * (ran->uniformRv() - 0.5) );
	setEventRate(cterm*oldEventRate);
	
	
	double LogPriorRatio = ran->lnExponentialPdf(poissonRatePrior, getEventRate()) - ran->lnExponentialPdf(poissonRatePrior, oldEventRate);
	double logProposalRatio = log(cterm);
	double logHR = LogPriorRatio + logProposalRatio;
	const bool acceptMove = acceptMetropolisHastings(logHR);
	
	//cout << "ER " << oldEventRate << "\t" << cterm*oldEventRate << endl;
	
	if (acceptMove == true){
		// continue
		acceptCount++;
		acceptLast = 1;
	}else{
		setEventRate(oldEventRate);
		rejectCount++;
		acceptLast = 0;
	}
	
	incrementGeneration();
	//cout << "Leaving UpdateEventRate" << endl;
}

 

void TraitModel::updateBetaMH(void){

	
#ifdef DEBUG_LH	
	double startLH = getCurrLnLTraits();
	vector<double> t1;
	vector<double> t2;
	
	for (int i = 0; i < treePtr->getNumberOfNodes(); i++)
		t1.push_back(treePtr->getNodeFromDownpassSeq(i)->getMeanBeta());
#endif		
	// Tempt stuff at start...
	
	
	//double startLL = getCurrLnLTraits();
	
 	int toUpdate = ran->sampleInteger(0, eventCollection.size());
	
	TraitBranchEvent* be = rootEvent;
	
	
	if (toUpdate > 0){
		std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
		for (int i = 1; i < toUpdate; i++){
			myIt++;		
		}
		
		be = (*myIt);
	}
	
	
	double oldRate = be->getBetaInit();
	double cterm = exp( _updateBetaScale * (ran->uniformRv() - 0.5) );
	be->setBetaInit(cterm * oldRate);
 	treePtr->setMeanBranchTraitRates();
	

	
	
	
#ifdef NO_DATA
	double PropLnLik = 0;
#else
	
	double PropLnLik = computeLikelihoodTraits();
	
#endif		
	
	double LogPriorRatio = ran->lnExponentialPdf(sttings->getBetaInitPrior(), be->getBetaInit());
	LogPriorRatio -= ran->lnExponentialPdf(sttings->getBetaInitPrior(), oldRate);
	
	double LogProposalRatio = log(cterm);
	
	double likeRatio = PropLnLik - getCurrLnLTraits();
	
	double logHR = likeRatio + LogPriorRatio + LogProposalRatio;
	
	const bool acceptMove = acceptMetropolisHastings(logHR);
	
//	cout << getGeneration() << "\tL1: " << startLH << "\tL2: " << getCurrLnLTraits() << endl;
	
	
	if (acceptMove == true){
		//cout << "accept: " << oldRate << "\t" << be->getBetaInit() << endl;
		setCurrLnLTraits(PropLnLik);
		acceptCount++;
		acceptLast = 1;
		
	}else{
		
		// revert to previous state
		_lastLH = PropLnLik;
		
#ifdef DEBUG_LH
		cout << "Before reset: " << computeLikelihoodTraits() << "\tAfter: ";
#endif		
		
		be->setBetaInit(oldRate);
 		treePtr->setMeanBranchTraitRates();
		acceptLast = 0;
		rejectCount++;	
		
#ifdef DEBUG_LH
		cout << computeLikelihoodTraits() << endl;
				
		double L1 = getCurrLnLTraits();
		double L2 = computeLikelihoodTraits();
		double delta = L1 - L2;
		
		for (int i = 0; i < treePtr->getNumberOfNodes(); i++)
			t2.push_back(treePtr->getNodeFromDownpassSeq(i)->getMeanBeta());		
		
		if ((delta < -0.01) || (delta > 0.01)){
			cout << "updateBetaMH:";
			cout << "problem at " << getGeneration() << "\tL1: " << getCurrLnLTraits();
			cout << "\tL2: " << computeLikelihoodTraits() << endl;
			exit(1);
	
		} 
#endif			
		
		
		
	}

	/*if (!acceptMove){
		cout << endl;
		cout << startLL << "\tCurr: " << getCurrLnLTraits() << "\tcalc: " << computeLikelihoodTraits() << endl;
 	}*/
	
	incrementGeneration();		
	
}


void TraitModel::updateBetaShiftMH(void){

	int toUpdate = ran->sampleInteger(0, eventCollection.size());
	
	TraitBranchEvent* be = rootEvent;
	
	
	if (toUpdate > 0){
		std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
		for (int i = 1; i < toUpdate; i++){
			myIt++;		
		}
		
		be = (*myIt);
	}
	
	double oldShift = be->getBetaShift();
	double newShift = oldShift + ran->normalRv((double)0.0, _updateBetaShiftScale);
 
	be->setBetaShift(newShift);
 	treePtr->setMeanBranchTraitRates();
	
#ifdef NO_DATA
	double PropLnLik = 0;
#else
	double PropLnLik = computeLikelihoodTraits();
	
#endif		
	
	// Normal prior on shift parameter:
	double LogPriorRatio = ran->lnNormalPdf((double)0.0, sttings->getBetaShiftPrior(), newShift);
	LogPriorRatio -= ran->lnNormalPdf((double)0.0, sttings->getBetaShiftPrior(), oldShift);

	
	double LogProposalRatio = 0.0;
	
	double likeRatio = PropLnLik - getCurrLnLTraits();
	
	double logHR = likeRatio + LogPriorRatio + LogProposalRatio;
	
	const bool acceptMove = acceptMetropolisHastings(logHR);
	
	if (acceptMove == true){
		
		setCurrLnLTraits(PropLnLik);
		acceptCount++;
		acceptLast = 1;
		
	}else{

		// revert to previous state
		be->setBetaShift(oldShift);
		treePtr->setMeanBranchTraitRates();
  		acceptLast = 0;
		rejectCount++;

		
#ifdef DEBUG_LH
		
		double L1 = getCurrLnLTraits();
		double L2 = computeLikelihoodTraits();
		double delta = L1 - L2;
		
		if ((delta < -0.01) || (delta > 0.01)){
			cout << "udpateBeatSHiftMH:";
			cout << "problem at " << getGeneration() << "\tL1: " << getCurrLnLTraits();
			cout << "\tL2: " << computeLikelihoodTraits() << endl;
			exit(1);
		} 
#endif			
		
	
	}
	
	
	
	incrementGeneration();	


}


void TraitModel::updateNodeStateMH(void){

	Node * xnode = treePtr->chooseInternalNodeAtRandom();
	
#ifdef NO_DATA
	
#else
	double oldTriadLogLik = computeTriadLikelihoodTraits(xnode);
	
#endif
	
	double oldstate = xnode->getTraitValue();
	double newstate = oldstate + ran->uniformRv((-1.0 * sttings->getUpdateNodeStateScale()), sttings->getUpdateNodeStateScale());
	xnode->setTraitValue(newstate);
	
#ifdef NO_DATA
	double PropLnLik = 0.0;
#else
	// Compute triad likelihood
	double newTriadLoglik = computeTriadLikelihoodTraits(xnode);
	double PropLnLik = getCurrLnLTraits() - oldTriadLogLik + newTriadLoglik;
	
#endif		
	
	// set to zero for now...flat prior, unifor (see below).
	double LogPriorRatio = 0.0;  
	
	
	double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0
	
	double likeRatio = PropLnLik - getCurrLnLTraits();
	
	
	double logHR = LogPriorRatio + logProposalRatio + likeRatio;
	bool acceptMove = acceptMetropolisHastings(logHR);		
	
	// Here we do prior calculation to avoid computing infinity...
	if (newstate > sttings->getTraitPriorMax() || newstate < sttings->getTraitPriorMin()){
		acceptMove = false;
	}
	if (acceptMove == true){
 		//continue
		setCurrLnLTraits(PropLnLik);		
		acceptCount++;
		acceptLast = 1;
		
	}else{
 		xnode->setTraitValue(oldstate);
		
		rejectCount++;
		acceptLast = 0;
	}
	
	incrementGeneration();	
	
	
}


void TraitModel::updateNodeStateMH(Node * xnode){
 
#ifdef NO_DATA
	
#else
	double oldTriadLogLik = computeTriadLikelihoodTraits(xnode);
	
#endif
	
	double oldstate = xnode->getTraitValue();
	double newstate = oldstate + ran->uniformRv((-1.0 * sttings->getUpdateNodeStateScale()), sttings->getUpdateNodeStateScale());
	xnode->setTraitValue(newstate);
	
#ifdef NO_DATA
	double PropLnLik = 0.0;
#else
	// Compute triad likelihood
	double newTriadLoglik = computeTriadLikelihoodTraits(xnode);
	double PropLnLik = getCurrLnLTraits() - oldTriadLogLik + newTriadLoglik;
	
	
	
#endif		
	
	// set to zero for now...flat prior, unifor (see below).
	double LogPriorRatio = 0.0;  
	
	
	double logProposalRatio = 0.0; // proposal ratio for uniform = 1.0
	
	double likeRatio = PropLnLik - getCurrLnLTraits();
	
	
	double logHR = LogPriorRatio + logProposalRatio + likeRatio;
	bool acceptMove = acceptMetropolisHastings(logHR);		
	
	// Here we do prior calculation to avoid computing infinity...
	if (newstate > sttings->getTraitPriorMax() || newstate < sttings->getTraitPriorMin()){
		acceptMove = false;
	}
	if (acceptMove == true){
		//continue
		setCurrLnLTraits(PropLnLik);		
		acceptCount++;
		acceptLast = 1;
		
	}else{
 		xnode->setTraitValue(oldstate);
		
		rejectCount++;
		acceptLast = 0;
	}
	
	incrementGeneration();	
	
	
	


}


void TraitModel::updateDownstreamNodeStatesMH(Node * xnode){

	// Get list of internal node descendants from node * x
	// update each (or some fraction thereof).
	
	treePtr->setTempInternalNodeArray(xnode);
	for (int i = 0; i < 100; i++){
		updateNodeStateMH(treePtr->getRandomNodeFromTempArray());
	}
	
	treePtr->clearTempNodeArray();


}


double TraitModel::computeLikelihoodTraits(void){

	double LnL = 0.0;
	
	//Node * tmpnode = treePtr->getRoot()->getLfDesc();
	
#ifdef NO_DATA
	LnL = 0.0;
#else
	int numNodes = treePtr->getNumberOfNodes();
	
	// iterate over non-root nodes and compute LnL
	
	for (int i = 0; i < numNodes; i++){
		Node * xnode = treePtr->getNodeFromDownpassSeq(i);
		if ( (xnode != treePtr->getRoot()) && (xnode->getCanHoldEvent() == true) ){
			
			
			double var = xnode->getBrlen() * xnode->getMeanBeta();
			
			// change in phenotype:
			double delta = xnode->getTraitValue() - xnode->getAnc()->getTraitValue();
 
			LnL += ran->lnNormalPdf(0, var, delta);
			
			//cout << xnode << "dz: " << delta << "\tT: " << xnode->getBrlen() << "\tRate: " << xnode->getMeanBeta();
			//cout << "\tLf: " << ran->lnNormalPdf(0, var, delta) << endl;
			
			/*if (xnode == tmpnode){
				cout << tmpnode->getTraitBranchHistory()->getAncestralNodeEvent()->getBetaInit();
				cout << "\tDelta: " << delta << "\tvar: " << var << "\tLL: " << ran->lnNormalPdf(0, var, delta);
				cout << "\tBeta: " << xnode->getMeanBeta()	<< endl;
			}*/
		}
		
	}
	
#endif		
	
	return LnL;

}

double TraitModel::computeTriadLikelihoodTraits(Node * x){
	
#ifdef DEBUG
	cout << "Enter computeTriadLikelihood: Node : " << x << endl;
#endif	
	
	double logL = 0.0;
	
	// Can only use this likelihood if node contributes to
	// likelihood of observed data
	
	if (x->getCanHoldEvent() == true){
		
		// computation for left descendant branch:
		
		if (x->getLfDesc()->getCanHoldEvent() == true){
			double delta = x->getLfDesc()->getTraitValue() - x->getTraitValue();
			logL += ran->lnNormalPdf(0, (x->getLfDesc()->getBrlen() * x->getLfDesc()->getMeanBeta()), delta);		
		}	
		
		
		if (x->getRtDesc()->getCanHoldEvent() == true){
			// computation for right descendant branch
			double delta = x->getRtDesc()->getTraitValue() - x->getTraitValue();
			logL += ran->lnNormalPdf(0, (x->getRtDesc()->getBrlen() * x->getRtDesc()->getMeanBeta()), delta);		
			
			
		}
		
		
		// computation for ancestral branch (unless == root)
		
		if (x != treePtr->getRoot()){
			
			double delta = x->getTraitValue() - x->getAnc()->getTraitValue();
			logL += ran->lnNormalPdf(0, (x->getBrlen() * x->getMeanBeta()), delta);	
			
		}	
		
		
		
	}
	
#ifdef DEBUG
	cout << "Leaving computeTriadLikelihood: Node : " << x << endl;
#endif	
	
	return logL;

}




double TraitModel::computeLogPrior(void){
	

	
	double logPrior = 0.0;
	logPrior += ran->lnExponentialPdf(sttings->getBetaInitPrior(), rootEvent->getBetaInit());
	logPrior += ran->lnNormalPdf((double)0.0, sttings->getBetaShiftPrior(), rootEvent->getBetaShift());	
	for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin(); i != eventCollection.end(); i++){
		logPrior += ran->lnExponentialPdf(sttings->getBetaInitPrior(), (*i)->getBetaInit());
		logPrior += ran->lnNormalPdf((double)0.0, sttings->getBetaShiftPrior(), (*i)->getBetaShift());
	}
 
	// and prior on number of events:
	
	logPrior += ran->lnPoissonProb((1/(double)sttings->getTargetNumberOfEvents()), eventCollection.size());
	
	return logPrior;
	
}





bool TraitModel::acceptMetropolisHastings(const double lnR) {
	const double r = safeExponentiation(TraitModel::mhColdness*lnR);
	return (ran->uniformRv() < r);
}



void TraitModel::initializeBranchHistories(Node* x){
	//cout << x << endl;
	x->getTraitBranchHistory()->setNodeEvent(rootEvent);
	
	if (x->getAnc() != NULL){
		x->getTraitBranchHistory()->setAncestralNodeEvent(rootEvent);
	}
	
	if (x->getLfDesc() != NULL){
		initializeBranchHistories(x->getLfDesc());
	}
	if (x->getRtDesc() != NULL){
		initializeBranchHistories(x->getRtDesc());
	}
	
	
}



void TraitModel::printStartAndEndEventStatesForBranch(Node* x){
	
	if (x != treePtr->getRoot()){
		cout << "Node: " << x << "\tAnc: " << x->getTraitBranchHistory()->getAncestralNodeEvent();
		cout << "\tevent: " << x->getTraitBranchHistory()->getNodeEvent() << endl;
	}
	
	if (x->getLfDesc() != NULL){
		printStartAndEndEventStatesForBranch(x->getLfDesc());
	}
	if (x->getRtDesc() != NULL){
		printStartAndEndEventStatesForBranch(x->getRtDesc());
	}
}




/*
 If this works correctly, this will take care of the following:
 1. if a new event is created or added to tree,
 this will forward set all branch histories from the insertion point
 
 2. If an event is deleted, you find the next event rootwards,
 and call forwardSetBranchHistories from that point. It will replace
 settings due to the deleted node with the next rootwards node.
 
 */

void TraitModel::forwardSetBranchHistories(TraitBranchEvent* x){
	// If there is another event occurring more recent (closer to tips)
	//	do nothing. Even just sits in branchHistory but doesn't affect
	//	state of any other nodes.
	
	
	// this seems circular, but what else to do?
	//	given an event (which references the node defining the branch on which event occurs)
	//	 you get the corresponding branch history and the last event
	//	 since the events will have been inserted in the correct order.
	
	Node* myNode = x->getEventNode();
	//cout << "Node: " << myNode << endl;
	
	//cout << endl << endl;
	//cout << "event in forwardSet: " << x << endl;
	
	//printEventData();
	
	
	
	
	if (x == rootEvent){
		forwardSetHistoriesRecursive(myNode->getLfDesc());
		forwardSetHistoriesRecursive(myNode->getRtDesc());	
		
	}else if (x == myNode->getTraitBranchHistory()->getLastEvent()){
		// If TRUE, x is the most tip-wise event on branch.
		myNode->getTraitBranchHistory()->setNodeEvent(x);
		
		// if myNode is not a tip:
		if (myNode->getLfDesc() != NULL && myNode->getRtDesc() != NULL){
			forwardSetHistoriesRecursive(myNode->getLfDesc());
			forwardSetHistoriesRecursive(myNode->getRtDesc());
		}
		// else: node is a tip : do nothing.
		
		
	}
	//else: there is another more tipwise event on same branch; do nothing 
	
	
}


void TraitModel::forwardSetHistoriesRecursive(Node* p){
	
	// Get event that characterizes parent node
	TraitBranchEvent* lastEvent = p->getAnc()->getTraitBranchHistory()->getNodeEvent();
	// set the ancestor equal to the event state of parent node:
	p->getTraitBranchHistory()->setAncestralNodeEvent(lastEvent);
	
	// if no events on the branch, go down to descendants and do same thing
	//	otherwise, process terminates (because it hits another event on branch
	if (p->getTraitBranchHistory()->getNumberOfBranchEvents() == 0){
		p->getTraitBranchHistory()->setNodeEvent(lastEvent);
		if (p->getLfDesc() != NULL){
			forwardSetHistoriesRecursive(p->getLfDesc());		
		}
		if (p->getRtDesc() != NULL){
			forwardSetHistoriesRecursive(p->getRtDesc());
		}
	}
	
}




void TraitModel::printBranchHistories(Node * x){
	
	if (x != treePtr->getRoot()){
		cout << "Node: " << x;
		cout << "\t#Events: " << x->getTraitBranchHistory()->getNumberOfBranchEvents() << "\tStart: ";
		cout << x->getTraitBranchHistory()->getAncestralNodeEvent() << "\tEnd: ";
		cout <<	x->getTraitBranchHistory()->getNodeEvent() << endl;
		
	}
	if (x->getLfDesc() != NULL)
		printBranchHistories(x->getLfDesc());
	if (x->getRtDesc() != NULL)
		printBranchHistories(x->getRtDesc());
	
	
	
}



double	TraitModel::getMHacceptanceRate(void){
	
	double arate = (double)acceptCount / ((double)acceptCount + (double)rejectCount);
	acceptCount = 0;
	rejectCount = 0;
	
	return arate;
	
}


TraitBranchEvent* TraitModel::getEventByIndex(int x){
	
	int ctr = 0;
	std::set<TraitBranchEvent*>::iterator myIt = eventCollection.begin();
	for (int i = 0; i <= x; i++){
		myIt++;
	}
	
	return (*myIt);
}



/*
 Model::countTimeVaryingRatePartitions
 
 -counts number of time-varying rate partitions
 
 */
int	TraitModel::countTimeVaryingRatePartitions(void){
	
	int count = 0;
	count += (int)rootEvent->getIsEventTimeVariable();
	for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin(); i != eventCollection.end(); i++){
		count += (int)(*i)->getIsEventTimeVariable();
	}
	return count;
}


/*
 Write event data to file for all events "on" tree
 at a given point in the MCMC chain
 
 
 */

void TraitModel::getEventDataString(stringstream &ss){
 
	
	ss << "GEN=" << getGeneration() << "[";
	
	
	TraitBranchEvent * be = rootEvent;
	Node * xl = treePtr->getRoot()->getRandomLeftTipNode();
	Node * xr = treePtr->getRoot()->getRandomRightTipNode();
	ss << "{taxon1=" << xl->getName() << ",taxon2=" << xr->getName() << ",time=" << be->getAbsoluteTime() << ",";
	ss << "beta1=" << be->getBetaInit() << ",beta2=" << be->getBetaShift() << "}";
	
	
	
	
	if (eventCollection.size() > 0){
		for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin(); i != eventCollection.end(); i++){
			be = (*i);
			if (be->getEventNode()->getLfDesc() == NULL){
				//cout << be->getEventNode()->getName() << endl;
				//cout << "xleft: " << xl->getName() << "\txrt: " << xr->getName() << endl;
				ss << "{taxon1=" << be->getEventNode()->getName() << ",taxon2=NA" << ",time=" << be->getAbsoluteTime() << ",";
				ss << "beta1=" << be->getBetaInit() << ",beta2=" << be->getBetaShift() << "}";			}else{
				Node * xl = be->getEventNode()->getRandomLeftTipNode();
				Node * xr = be->getEventNode()->getRandomRightTipNode();
				
				//cout << "xleft: " << xl->getName() << "\txrt: " << xr->getName() << endl;
				ss << "{taxon1=" << xl->getName() << ",taxon2=" << xr->getName() << ",time=" << be->getAbsoluteTime() << ",";
				ss << "beta1=" << be->getBetaInit() << ",beta2=" << be->getBetaShift() << "}";			}
			
			
			
		}	
		
	}
	
	
	ss << "];";
	
}


bool TraitModel::isEventConfigurationValid(TraitBranchEvent * be){
	//cout << "enter isEventConfigValid" << endl;
	bool isValidConfig = false;
	
	if (be->getEventNode() == treePtr->getRoot()){
		Node * rt = treePtr->getRoot()->getRtDesc();
		Node * lf = treePtr->getRoot()->getLfDesc();
		if (rt->getTraitBranchHistory()->getNumberOfBranchEvents() > 0 && lf->getTraitBranchHistory()->getNumberOfBranchEvents() > 0){
			// events on both descendants of root. This fails.
			isValidConfig = false;
		}else{
			isValidConfig = true;
		}
		
	}else{
		int badsum = 0;
		
		Node * anc = be->getEventNode()->getAnc();
		Node * lf = anc->getLfDesc();
		Node * rt = anc->getRtDesc();
		
		//cout << "a: " << anc << "\tb: " << lf << "\tc: " << rt << endl;
		
		// test ancestor for events on branch:
		
		if (anc == treePtr->getRoot()){
			badsum++;
		}else if (anc->getTraitBranchHistory()->getNumberOfBranchEvents() > 0){
			badsum++;
		}else{
			// nothing;
		}
		
		// test lf desc:
		if (lf->getTraitBranchHistory()->getNumberOfBranchEvents() > 0){
			badsum++;
		}
		
		// test rt desc
		if (rt->getTraitBranchHistory()->getNumberOfBranchEvents() > 0){
			badsum++;
		}
		
		if (badsum == 3){
			isValidConfig = false;
		}else if (badsum < 3){
			isValidConfig = true;
		}else{
			cout << "problem in Model::isEventConfigurationValid" << endl;
			exit(1);
		}
		
		
	}
	
	
	//cout << "leaving isEventConfigValid. Value: " << isValidConfig << endl;
	return isValidConfig;
}



void TraitModel::printEventData(void){
	
	TraitBranchEvent * be = rootEvent;
	cout << "RtBt: " << be->getBetaInit() << "\tSf: " << be->getBetaShift() << "\t";
	int ctr = 0;
	for (std::set<TraitBranchEvent*>::iterator i = eventCollection.begin(); i != eventCollection.end(); i++){
		be = (*i);
		cout << ctr++ << "\tBt: " << be->getBetaInit() << "\tSt: " << be->getBetaShift() << "\tMap: " << be->getMapTime() << endl;;
	
	}
	cout << endl;
}

/*
void TraitModel::initializeTraitParamsForNodes(void){

	for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++){
		
	
	}


}



*/


void TraitModel::setMinMaxTraitPriors(void){
	
	int nnodes = treePtr->getNumberOfNodes();
	vector<double> tvec;
	for (int i = 0; i < nnodes; i++){
		Node * xnode = treePtr->getNodeFromDownpassSeq(i);
		if (xnode->getTraitValue() != NULL){
			tvec.push_back(xnode->getTraitValue());
		}
	}
	
	sort(tvec.begin(), tvec.end());
	
	// cout << "Min: " << tvec[0] << "\tMax: " << tvec[(tvec.size() - 1)] << endl;
	
	// Default here will be to use observed range +/- 20% 
	double rg = tvec[(tvec.size() - 1)] - tvec[0];
	double minprior = tvec[0] - (0.2 * rg);
	double maxprior = tvec[(tvec.size() - 1)] + (0.2 * rg);
	
	cout << "Min and max phenotype limits set using observed data: " << endl;
	cout << "\t\tMin: " << minprior << "\tMax: " << maxprior << endl;
	sttings->setTraitPriorMin(minprior);
	sttings->setTraitPriorMax(maxprior);
	
}



















