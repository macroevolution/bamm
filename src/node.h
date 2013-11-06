/*
 *  node.h
 *  proj7
 *
 *  Created by Dan Rabosky on 9/16/11.
  *
 */

#ifndef Node_H
#define Node_H

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


class Tree{

private:
	Node*			root;
	set<Node*>		nodes;
	vector<Node*>	downPassSeq;
	
	// Internal node set:: for choosing random node to update state
	set<Node*>		internalNodeSet;
	
	MbRandom*		ranPtr;
	double			_startTime;
	double			_tmax;
	bool			_isExtant;
	int				_ntaxa;
	void			recursivelyAddNodesToSet(Node * p);
	void			rebuildTreeNodeSet(void);
	double			_age;	// time to root node, from present
	
	void			setIsLivingTipStatus(void);
	void			getDownPassSeq(void);
	void			passDown(Node * p);
	void			setTipStatus(void);
	
	double			treelength; // this is the total treelength: also used in mapping events to tree
	
	// this is a pointer to an object that stores 
	//	all events that happened.
	eventSet*		treeEvents; // What is this doing???
	
	// New data for mapping events to nodes:
	set<Node*>		mappableNodes; // Nodes subtended by a branch that can hold an event
	double			_totalMapLength;
	
	set<Node*>		_tempNodeSet;

	
public:

					Tree(string fname, MbRandom * rnptr);
	
					Tree(void);
					~Tree(void);
	// this function initializes eventHistory for tree, 
	//	a collection of pointers to events...
	
	double			getTreeLength(void)		{ return treelength;	}
	Node*			mapEventToTree(double x);	
	
	
	// this function will map branches such that each branch
	//has a unique segment of the line defined on (O, treelength). 
	
	void			setTreeMap(Node* p); 
	void			printNodeMap(void);
	void			setTreeMap(int ndesc); // set tree map for some number of descendant nodes...
 	double			getAbsoluteTimeFromMapTime(double x);
	
	
	Node*			getDownPassNode(int i)	{ return downPassSeq[i]; }
	int				getNumberOfNodes(void)	{ return (int)downPassSeq.size(); }	
	Node*			getNodeFromDownpassSeq(int i)	{ return downPassSeq[i]; }
	
	// count number of descendant nodes from a given node
	int				getDescNodeCount(Node* p);
	int				getDescTipCount(Node* p); // get number of tips from given node
	int				countDescendantsWithValidTraitData(Node * p);
	
	// If all nodes are valid eg for speciation extinction model, 
	//	just call this, which sets all node status to TRUE
	void			setAllNodesCanHoldEvent(void);
	void			setCanNodeHoldEventByDescCount(int x);
	
					/* disposable stuff below? */ 
	void			setStartTime(double x)	{ _startTime = x; }
	void			setTmax(double x)	{ _tmax = x; }
	double			getTmax(void)	{ return _tmax; }
	//Node*			chooseRandomNode(set<Node*> nodeList);
	
	Node*			getRoot(void)		{ return root; }
	
	string			getNewick(void);
	void			writeTree(Node * p, stringstream &ss);
	int				getNumberTips(void);
	int				getNumberExtantTips(void);
	void			pruneExtinctTaxa(void); // removed
	void			fixExtinct(Node * p);
	void			setExtantStatus(bool x)	{ _isExtant = x;	}
	void			setExtantStatus(void);
	bool			getExtantStatus(void)		{ return _isExtant; }
	void			writeNodeData(void);
	void			setBranchLengths(void);
	void			deleteExtinctNodes(void);
	void			buildTreeFromNewickString(string ts);
	void			setTaxonCountFromNewickString(string ts);
	//int			getNtaxa(void)				{ return _ntaxa; } /*bad*/
	bool			isValidChar(char x);
	
	void			setNodeTimes(Node * p);	
	void			setBranchingTimes(Node * p);
	
	double			getAge(void) {  return _age; }
	void			setAge(void);
	vector<double>  getBranchingTimes(void);
	//double			computeGammaStat(void);
	//void			pruneRandomTaxa(int prunecount);	
 
	void			writeMeanBranchTraitRateTree(Node * p, stringstream &ss);
	void			setMeanBranchTraitRates(void);
	
	// Functions for phenotypic evolution:
	void			getPhenotypes(string fname);
	void			getPhenotypesMissingLatent(string fname);
	
	void			printTraitValues(void);
	void			initializeTraitValues(void);
	void			recursiveSetTraitValues(Node * x, double mn, double mx);
	Node*			chooseInternalNodeAtRandom(void);

	void			generateTraitsAllNodesBM(Node* xnode, double varx);
	void			generateTraitsAllNodesFromEventBeta(Node* xnode);
	void			printTraitRange(void);
	double			getTraitMinTip(void);
	double			getTraitMaxTip(void);
	
	// speciation-extinction output
	void			setMeanBranchSpeciation(void);
	void			setMeanBranchExtinction(void);

	void			echoMeanBranchRates(void);

	void			writeMeanBranchSpeciationTree(Node * p, stringstream &ss);
	void			writeBranchSpeciationRatesToFile(string fname, bool append);
	void			writeBranchExtinctionRatesToFile(string fname, bool append);
	
 
	void			writeMeanBranchExtinctionTree(Node * p, stringstream &ss);
	void			writeNodeSpeciationTree(Node * p, stringstream &ss);
	
	void			writeMeanBranchNetDivRateTree(Node * p, stringstream &ss);
	void			writeBranchPhenotypes(Node * p, stringstream & ss);
	
	// speciation-extinction initialization:
	void			initializeSpeciationExtinctionModel(string fname); // file for species-specific values.
	void			initializeSpeciationExtinctionModel(double sampFraction);
	void			printInitialSpeciationExtinctionRates(void);
	
	// New fxns for mapping events to nodes:
	
	double			getTotalMapLength(void)		{ return _totalMapLength;	}
	void			setTotalMapLength(double x)	{ _totalMapLength = x;		}	
	void			setCanNodeBeMapped(int ndesc);  

	
	void			loadPreviousNodeStates(Tree * ostree);
	
	
	// Functions for random access of nodes from temporary nodeset array
	
	void			setTempInternalNodeArray(Node * p);
	void			tempNodeSetPassDown(Node * p);
	void			clearTempNodeArray(void);
	Node*			getRandomNodeFromTempArray(void);
	
	
	void			printNodeBranchRates(void);	
	
	void			computeMeanTraitRatesByNode(Node * x);
	
	//set<Node*>		getInternalDescendantSet(Node * p);	
	//Node*			passDownReturnNode(Node * p, set<Node*> * nodeAr);
	Node *			getNodeMRCA(string A, string B);
	void			passUpFillTempNodeArray(Node * x);
	Node *			getNodeByName(string A);
 
	void			printNodeTraitRates(void);

	void			printCanHoldEventByNode(void);	
	

	void			echoMeanBranchTraitRates(void);
	
};

 

#endif







