

#ifndef MODEL_H
#define MODEL_H

#include <set>
#include <vector>
 
#include "event.h"

using namespace std;

//Forward declarations
class Tree;
class Node;
class MbRandom;
class Settings;

double safeExponentiation(double x);
double proportionalShrink(double x, double scale);

class Model
{
public:

	static double	mhColdness;	
	
	Model(MbRandom* ranptr, Tree * tp, Settings* sp);
	~Model(void);
	
	
	// Full likelihood will be lnLikTraits + lnLikBranches
	void			setCurrLnLTraits(double x)		{ lnLikTraits = x;		}
	double			getCurrLnLTraits(void)			{ return lnLikTraits;	}
	void			setCurrLnLBranches(double x)	{ lnLikBranches = x;	}
	double			getCurrLnLBranches(void)		{ return lnLikBranches; }
	
	// Full likelihood function
 
	
	// Likelihood functions for branches data

	double			computeLikelihoodBranches(void);
	double			computeLikelihoodBranchesByInterval(void);
	
	
//	double			computeLikelihoodBranchesByBranchMeanRates(void);
//	Decremented: 9.20.2012  * see previous version
	
	// likelihood functions for trait data
	//double			computeLikelihoodTraits(void);
	//double			computeTriadLikelihoodTraits(Node * x);
	
	double			computeLogPrior(void);
	
/* deprecated */ // double			verboseComputeLikelihood(void);
	
	
	double			getEventRate(void)		{ return eventLambda; }
	void			setEventRate(double x)	{ eventLambda = x;    }
	
	double			proportionalShrink(double x, double scale);	
	bool			acceptMetropolisHastings(const double lnR);
	void			incrementGeneration(void)	{	gen++;	}
	int				getGeneration(void)			{ return gen;  }
	void			resetGeneration(void)		{ gen = 0;	} // to be used after TraitPreBurnin
	

// Stuff for event handling:
	void			addEventToTree(void); 
	void			addEventToTree(double x); // add to specific map position...
	void			deleteRandomEventFromTree(void);
	void			deleteEventFromTree(BranchEvent * be);
	
	void			printEvents(void); 

	int				getNumberOfEvents(void)		{ return eventCollection.size(); }
		
	
	BranchEvent*	getRootEvent(void)			{ return rootEvent; }
	

	// These functions take a branch event
	//	and recursively update branch histories for all nodes
	//	going towards the tips
	void			forwardSetBranchHistories(BranchEvent* x);
	void			forwardSetHistoriesRecursive(Node* p);
 
	int				countEventsInBranchHistory(Node* p);
	
	
//	BranchEvent*	getLastEvent(Node* p);
//	BranchEvent*	getLastEvent(BranchEvent* x); 
	
	// initialize all branch histories to the root node.
	void			initializeBranchHistories(Node* x);
	//void			printBranchHistories(Node* x);	
	void			printStartAndEndEventStatesForBranch(Node* x);
	
	
	void			eventLocalMove(BranchEvent* x); // move specific event
	void			eventGlobalMove(BranchEvent* x); // move specific event
	void			eventLocalMove(void); // move random event
	void			eventGlobalMove(void); // move random event
	void			revertEventToPreviousPosition(void);
	
	
	
	
	// Return random event, or NULL if no events on tree (other than root)
	BranchEvent*	chooseEventAtRandom(void); 
	
	// check whether branch histories are set correctly recursively
//	void			checkAreBranchHistoriesSetCorrectly(Node * p);
//	void			checkEventSetAgainstBranchHistories(Node * p);

	// MCMC:
	void			changeNumberOfEventsMH(void);	// propose addition or deletion; accept/reject move.
	void			moveEventMH(void);
	void			revertMovedEventToPrevious(void);
	
	/*	*****************	*/ 
	//void			updateSpeciationRateMH(void);
	//void			updateBetaMH(void);
	// lambda/mu related stuff:
	void			updateLambdaInitMH(void);
	void			updateLambdaShiftMH(void);

	void			updateMuInitMH(void);
	void			updateMuShiftMH(void);
	
	void			updateNodeStateMH(void);
	void			updateNodeStateMH(Node * xnode);
	void			updateDownstreamNodeStatesMH(Node * xnode);
	void			updateEventRateMH(void);
	void			updateTimeVariablePartitionsMH(void);

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
	
	BranchEvent*	getEventByIndex(int x);			
	
	void			printExtinctionParams(void);
	int				countTimeVaryingRatePartitions(void);
	
	// Generate string with event data:
	void			getEventDataString(stringstream &ss);
	
	bool			isEventConfigurationValid(BranchEvent * be);
	
	void			initializeModelFromEventDataFile(void);
	
	void			debugLHcalculation(void);
	
private:
	
//	parameters of the model:
 
	// deprecated in favor of:
	//double			lnLik;
	double			lnLikTraits;
	double			lnLikBranches;
	
	// new parameters: March 23 2012
	double			_updateLambdaInitScale;
	double			_updateLambdaShiftScale;
	
	double			_updateMuInitScale;
	double			_updateMuShiftScale;
	
	// This parameter holds the density of the new 
	// parameters proposed during jump moves.
	// If the parameters are sampled from the prior,
	//		these should exactly cancel.
	
	double			_logQratioJump;
	
	// Root event parameters:
	
	double			_lambdaInit0;
	double			_lambdaShift0;
	
	double			_muInit0;
	double			_muShift0;
	
	
	// Parameters for MCMC proposals
	double				_scale; // scale for moving event
	double				_updateEventRateScale;
	double				_targetNumber;
	double				_localGlobalMoveRatio;
	
	double				eventLambda; // Poisson rate
	//double				GainLossRatio;
	//double				localGlobalMoveRatio; // frequency of local vs global event moves
	
	// Priors
	//double				betaPrior; // exponential rate
	double				poissonRatePrior; // exponential  /* keep this in Model */
	double				_lambdaInitPrior;
	double				_lambdaShiftPrior;
	
	double				_muInitPrior;
	double				_muShiftPrior;
	
	// Other private variables
	
	int					gen;
	MbRandom*			ran;
	Tree*				treePtr;
	Settings*			sttings;
 
	int					acceptCount;
	int					rejectCount;
	
	
	set<BranchEvent*>	eventCollection; // does NOT contain root event
	BranchEvent*		rootEvent; //branch event at root node; can't be modified
	double				lastDeletedEventMapTime; // map time of last deleted event
	
	double				_lastDeletedEventLambdaInit;
	double				_lastDeletedEventLambdaShift;
	
	double				_lastDeletedEventMuInit;
	double				_lastDeletedEventMuShift;
 
	
	
	// HEre are several variables that track the previous
	// state. At some point, these should have their own class
	
	// this is a pointer to the last event modified, whether 
	// it is moved, or has lambda updated, or whatever.
	BranchEvent*		lastEventModified;
	
	// General acceptreject flag:
	int					acceptLast; // true if last generation was accept; false otherwise
 
	double				_segLength; // for splitting branches
	
	
};


#endif





