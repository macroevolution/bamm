#ifndef TREE_H
#define TREE_H

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
class Node;


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
