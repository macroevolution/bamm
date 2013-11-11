/*
 *  Node.h
 *  proj7
 *
 *  Created by Dan Rabosky on 9/16/11.
  *
 */

#ifndef NODE_H
#define NODE_H

#include <string>
#include <set>
#include <vector>

using namespace	std;

class MbRandom;
class branchEvent;
class eventSet;
class Phenotype;
	
 
class BranchHistory;
class TraitBranchHistory;


class Node{

private:
	Node * _lfDesc;
	Node * _rtDesc;
	Node * _anc;
	string _name;
	
	string _cladeName;
	
	int _index;
	double _time;
	double _brlen; //
	double _branchTime;
	
	int _tipDescCount;
	bool _isExtant;
	bool _isTip;
	bool _isConstant;
	bool _isLivingTip;
	
	// specific for mapping:
	double _mapStart;
	double _mapEnd;
	
	BranchHistory*	history;
	TraitBranchHistory* _traitHistory;
	
	
	// For phenotypes:
	//Phenotype*		pheno;
	double	_trait; // trait value
	double	_meanBeta; // mean phenotypic rate
	double	_nodeBeta; // exact value at node.
	bool	_isTraitFixed; // is trait value a free parameter?
	
	
	
	//specific stuff for compound poisson rate model
	double _meanSpeciationRate;
	double _meanExtinctionRate;
	
	// Node rates for time-varying models:
	double _nodeLambda;
	double _nodeMu;
	
	
	
	double _di; // initial value of speciation probability (at node)
	double _ei; // initial value of extinction probability (at node)
	double _etip; // initial value (sampling fraction) at tip descended from node
	
	// value of log-likelihood at end of branch (after combining w speciation event):
	//		Zero if node is terminal. 
	double _nodeLikelihood; 

	// flag for whether node can or cannot define branch that can hold event:
	bool	_canHoldEvent;
	
 
	
public:
	Node(void); 
	Node(int x);
	void setLfDesc(Node * x)	{ _lfDesc = x; }
	void setRtDesc(Node * x)	{ _rtDesc = x; }
	void nullifyLfDesc(void)	{ _lfDesc = NULL; }
	void nullifyRtDesc(void)	{ _rtDesc = NULL; }
	void nullifyAnc(void)		{ _anc = NULL; }
	void setAnc(Node * x)		{ _anc = x; }
	void setName(string x)		{ _name = x; }
	void setIndex(int x)		{ _index = x; }
	void setTime(double x)		{ _time = x; }
	void setBrlen(double x)		{ _brlen = x; }
	
	void setTipDescCount(int x)	{ _tipDescCount = x; }
	
	
	Node * getLfDesc(void)		{ return _lfDesc; }
	Node * getRtDesc(void)		{ return _rtDesc; }
	Node * getAnc(void)			{ return _anc; }
	string getName(void)		{ return _name; }
	int getIndex(void)			{ return _index; }
	double getTime(void)		{ return _time; }
	double getBrlen(void)		{ return _brlen; }
	int getTipDescCount(void)	{ return _tipDescCount; }
	int getDescCount(void);
	void setExtantStatus(bool x) { _isExtant = x;	}  
	bool getExtantStatus(void)	{ return _isExtant; }
	
	void setIsTip(bool x)		{ _isTip = x;		}
	bool getIsTip(void)			{ return _isTip;	}
	void setIsConstant(bool x)	{ _isConstant = x; }
	bool getIsConstant(void)	{ return _isConstant; }
	
	void setIsLivingTip(bool x)	{ _isLivingTip = x; }
	bool getIsLivingTip(void)	{ return _isLivingTip;  }

	// Get a random TIP node descended from lfdesc of a given node
	Node *	getRandomLeftTipNode(void);

	// Get random RIGHT TIP node
	Node *	getRandomRightTipNode(void);
	
	// Specific for treemap:
	void setMapStart(double x)	{ _mapStart = x; }
	double getMapStart(void)	{ return _mapStart; }
	void setMapEnd(double x)	{ _mapEnd = x; }
	double getMapEnd(void)		{ return _mapEnd; }
	
	//Need to includet this
	BranchHistory* getBranchHistory(void)	{ return history;  }
	TraitBranchHistory* getTraitBranchHistory(void)	{ return _traitHistory;	}
	
	
	void	setMeanSpeciationRate(double x)	{ _meanSpeciationRate = x; }
	double	getMeanSpeciationRate(void)		{  return _meanSpeciationRate;  }
	
	void	setMeanExtinctionRate(double x)	{ _meanExtinctionRate = x;		}
	double	getMeanExtinctionRate(void)		{ return _meanExtinctionRate;	}
	
	void	computeNodeBranchSpeciationParams(void);
	void	computeNodeBranchExtinctionParams(void);
	
	void	setNodeLambda(double x)			{ _nodeLambda = x;		} 
	double	getNodeLambda(void)				{ return _nodeLambda;	}
	
	void	setNodeMu(double x)				{ _nodeMu = x;			}
	double	getNodeMu(void)					{ return _nodeMu;		}
	
	void	setNodeBeta(double x)			{ _nodeBeta = x;		}
	double	getNodeBeta(void)				{ return _nodeBeta;		}
	

	/*********/
	
	// phenotypic stuff:
	void	setTraitValue(double x)			{ _trait = x;		}
	double	getTraitValue(void)				{ return _trait;	}
	void	setMeanBeta(double x)			{ _meanBeta = x;	}
	double	getMeanBeta(void)				{ return _meanBeta; }
	void	setIsTraitFixed(bool x)			{ _isTraitFixed = x;}
	bool	getIsTraitFixed(void)			{ return _isTraitFixed; }
	
	// Speciation-extinction calculations
	void	setEinit(double x)					{ _ei = x;		}
	double	getEinit(void)						{ return _ei;	}
	void	setDinit(double x)					{ _di = x;		}
	double	getDinit(void)						{ return _di;	}
	void	setEtip(double x)					{ _etip = x;	}
	double	getEtip(void)						{ return _etip; }

	void	setNodeLikelihood(double x)			{ _nodeLikelihood = x;		}
	double	getNodeLikelihood(void)				{ return _nodeLikelihood;	}

	bool	getCanHoldEvent(void)			{ return _canHoldEvent; };
	void	setCanHoldEvent(bool x)			{ _canHoldEvent = x;	};

	// branching times
	void	setBranchTime(double x)			{ _branchTime = x;		}
	double	getBranchTime(void)				{ return _branchTime;	}
	
	void	setCladeName(string x)			{ _cladeName = x;		}
	string	getCladeName(void)				{ return _cladeName;	}

	//double	computeSpeciationRateInterval(double tstart, double tstop);
	double	computeSpeciationRateIntervalRelativeTime(double tstart, double tstop);
	double	computeSpeciationRateIntervalAbsoluteTime(double tstart, double tstop);
	double	computeExtinctionRateIntervalRelativeTime(double tstart, double tstop);
 
	
	double	getPointExtinction(double branchtime);
	
	
};


#endif
