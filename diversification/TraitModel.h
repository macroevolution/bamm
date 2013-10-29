/*
 *  TraitModel.h
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
 
 *
 */


#ifndef TRAITMODEL_H
#define TRAITMODEL_H

#include <set>
#include <vector>

#include "TraitEvent.h"

using namespace std;

//Forward declarations
class Tree;
class Node;
class MbRandom;
class Settings;
class TraitBranchEvent;

double safeExponentiation(double x);
double proportionalShrink(double x, double scale);

class TraitModel
{
public:
	
	static double	mhColdness;	
	
	TraitModel(MbRandom* ranptr, Tree * tp, Settings* sp);
	~TraitModel(void);
	
	// Full likelihood will be lnLikTraits + lnLikBranches
	void			setCurrLnLTraits(double x)		{ lnLikTraits = x;		}
	double			getCurrLnLTraits(void)			{ return lnLikTraits;	}
 
 
	double			computeLikelihoodTraits(void);
	double			computeTriadLikelihoodTraits(Node * x);
	
	double			computeLogPrior(void);
	
	double			getEventRate(void)		{ return eventLambda; }
	void			setEventRate(double x)	{ eventLambda = x;    }
	
	double			proportionalShrink(double x, double scale);	
	bool			acceptMetropolisHastings(const double lnR);
	void			incrementGeneration(void)	{	gen++;	}
	int				getGeneration(void)			{ return gen;  }
	void			resetGeneration(void)		{ gen = 0;	} // to be used after TraitPreBurnin
	
	
	// Stuff for event handling:
	void			addEventToTree(void); 
	void			addEventToTreeWithSetBeta(double beta, double bshift); 
	void			addEventToTree(double x); // add to specific map position...
	void			deleteRandomEventFromTree(void);
	void			deleteEventFromTree(TraitBranchEvent * be);
	
	void			printEvents(void); 
	
	int				getNumberOfEvents(void)		{ return eventCollection.size(); }
	
	TraitBranchEvent*	getRootEvent(void)			{ return rootEvent; }
	
	
	// These functions take a branch event
	//	and recursively update branch histories for all nodes
	//	going towards the tips
	void			forwardSetBranchHistories(TraitBranchEvent* x);
	void			forwardSetHistoriesRecursive(Node* p);
	
	int				countEventsInBranchHistory(Node* p);
 
	
	// initialize all branch histories to the root node.
	void			initializeBranchHistories(Node* x);
	//void			printBranchHistories(Node* x);	
	void			printStartAndEndEventStatesForBranch(Node* x);
	
	
	void			eventLocalMove(TraitBranchEvent* x); // move specific event
	void			eventGlobalMove(TraitBranchEvent* x); // move specific event
	void			eventLocalMove(void); // move random event
	void			eventGlobalMove(void); // move random event
	void			revertEventToPreviousPosition(void);
	
	
	
	
	// Return random event, or NULL if no events on tree (other than root)
	TraitBranchEvent*	chooseEventAtRandom(void); 
 
	
	// MCMC:
	void			changeNumberOfEventsMH(void);	// propose addition or deletion; accept/reject move.
	void			moveEventMH(void);
	void			revertMovedEventToPrevious(void);
	
	/*	*****************	*/ 
 
	// Trait evolution stuff:
	void			updateBetaMH(void);
	void			updateNodeStateMH(void);
	void			updateNodeStateMH(Node * xnode);
	void			updateBetaShiftMH(void);
	void			updateDownstreamNodeStatesMH(Node * xnode);
	void			updateEventRateMH(void);
	void			updateTimeVariablePartitionsMH(void);
	void			setMinMaxTraitPriors(void);
	
	// Probably Not necessary::
	//double			betaProposal(double traitrate);
	//double			getMeanTraitRate(void);
	
	/*	********************	*/
	
	
	
	void			printEventData(void);
	void			restoreLastDeletedEvent(void);
	
	// troubleshooting
	void			printBranchHistories(Node * x);
	
	Tree*			getTreePtr(void)		{ return treePtr; }
	
	// more output: acceptance rates
	double			getMHacceptanceRate(void);
	
	int				getAcceptLastUpdate(void)		{ return acceptLast;	}
	void			setAcceptLastUpdate(int x)		{ acceptLast = x;		}
	
	// 0 = last was rejected; 1 = accepted; -1 = not set.
	
	//double			getBetaPrior(void)			{ return betaPrior; }
	
	void			setPoissonRatePrior(double x)	{ poissonRatePrior	= x; }
	double			getPoissonRatePrior(void)		{ return poissonRatePrior;  }
	
	TraitBranchEvent*	getEventByIndex(int x);			
	
	int				countTimeVaryingRatePartitions(void);
	
	// Generate string with event data:
	void			getEventDataString(stringstream &ss);
	bool			isEventConfigurationValid(TraitBranchEvent * be);
	
	double			getLastLH(void)	{  return _lastLH;	}
	
private:
	
	//	parameters of the model:
	
	double				lnLikTraits;
	
	double				eventLambda; // Poisson rate
 
	double				poissonRatePrior; // exponential  /* keep this in Model */
	
	double				_updateBetaScale;
	double				_updateBetaShiftScale;
	double				_updateNodeStateScale;
	double				_scale;
	double				_updateEventRateScale;
	double				_localGlobalMoveRatio;
	double				_targetNumber;
	
	// Other private variables
	
	int					gen;
	MbRandom*			ran;
	Tree*				treePtr;
	Settings*			sttings;
	
	int					acceptCount;
	int					rejectCount;
	
	
	set<TraitBranchEvent*>	eventCollection; // does NOT contain root event
	TraitBranchEvent*		rootEvent; //branch event at root node; can't be modified
	double					lastDeletedEventMapTime; // map time of last deleted event
	
	double					_lastDeletedEventBetaInit;;
	double					_lastDeletedEventBetaShift;
 
	// HEre are several variables that track the previous
	// state. At some point, these should have their own class
	
	// this is a pointer to the last event modified, whether 
	// it is moved, or has lambda updated, or whatever.
	TraitBranchEvent*		lastEventModified;
	
	// General acceptreject flag:
	int						acceptLast; // true if last generation was accept; false otherwise
 
	// Ultimately, initializations should be handled in TraitModel and Model classes
	//	NOT in class TREE!!
	void			initializeTraitParamsForNodes(void);
	double			_lastLH;
	
	
};


#endif
