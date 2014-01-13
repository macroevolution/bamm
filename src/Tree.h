#ifndef TREE_H
#define TREE_H

#include <string>
#include <set>
#include <vector>

class MbRandom;
class branchEvent;
class eventSet;
class Phenotype;
class BranchHistory;
class TraitBranchHistory;
class Node;


// Branch lengths must be within this value to be considered the same
#define ULTRAMETRIC_ERROR 1e-4


class Tree
{

private:

    std::vector<double> terminalPathLengthsToRoot();
    void storeTerminalPathLengthsToRootRecurse
        (Node* node, std::vector<double>& pathLengths);
    bool allValuesAreTheSame(std::vector<double>& list);

    void assertTreeIsBifurcating();
    void assertTreeIsBifurcatingRecurse(Node* node);
    void assertBranchLengthsArePositive();
    void assertBranchLengthsArePositiveRecurse(Node* node);
    void assertTreeIsUltrametric();
    void assertTipsHaveUniqueNames();
    std::vector<std::string> terminalNames();
    void storeTerminalNamesRecurse(Node* node, std::vector<std::string>& names);

    Node* root;
    std::set<Node*> nodes;
    std::vector<Node*> downPassSeq;

    // Internal node set:: for choosing random node to update state
    std::set<Node*> internalNodeSet;

    MbRandom* ranPtr;
    double _startTime;
    double _tmax;
    bool _isExtant;
    int _ntaxa;
    void recursivelyAddNodesToSet(Node* p);
    void rebuildTreeNodeSet();
    double _age; // time to root node, from present

    void setIsLivingTipStatus();
    void getDownPassSeq();
    void passDown(Node* p);
    void setTipStatus();

    // This is the total treelength: also used in mapping events to tree
    double treelength;

    // This is a pointer to an object that stores
    // all events that happened.
    eventSet* treeEvents; // What is this doing???

    // New data for mapping events to nodes:
    // Nodes subtended by a branch that can hold an event
    std::set<Node*> mappableNodes;
    double _totalMapLength;

    std::set<Node*> _tempNodeSet;

public:

    Tree(std::string fname, MbRandom* rnptr);

    Tree();
    ~Tree();

    // This function initializes eventHistory for tree,
    // a collection of pointers to events...

    double getTreeLength();
    Node*  mapEventToTree(double x);

    // This function will map branches such that each branch
    // has a unique segment of the line defined on (O, treelength).

    void setTreeMap(Node* p);
    void printNodeMap();

    // Set tree map for some number of descendant nodes...
    void setTreeMap(int ndesc);

    double getAbsoluteTimeFromMapTime(double x);

    Node* getDownPassNode(int i);
    int   getNumberOfNodes();
    Node* getNodeFromDownpassSeq(int i);

    // Count number of descendant nodes from a given node
    int getDescNodeCount(Node* p);
    int getDescTipCount(Node* p); // get number of tips from given node
    int countDescendantsWithValidTraitData(Node* p);

    // If all nodes are valid eg for speciation extinction model,
    // just call this, which sets all node status to TRUE
    void setAllNodesCanHoldEvent();
    void setCanNodeHoldEventByDescCount(int x);

    /* disposable stuff below? */
    void   setStartTime(double x);
    void   setTmax(double x);
    double getTmax();

    Node* getRoot();

    std::string getNewick();
    void   writeTree(Node* p, std::stringstream& ss);

    int  getNumberTips();
    int  getNumberExtantTips();
    void pruneExtinctTaxa(); // removed
    void fixExtinct(Node* p);
    void setExtantStatus(bool x);
    void setExtantStatus();
    bool getExtantStatus();

    void writeNodeData();
    void setBranchLengths();
    void deleteExtinctNodes();
    void buildTreeFromNewickString(std::string ts);
    void setTaxonCountFromNewickString(std::string ts);
    bool isValidChar(char x);

    void setNodeTimes(Node* p);
    void setBranchingTimes(Node* p);

    double getAge();
    void   setAge();
    std::vector<double> getBranchingTimes();

    void writeMeanBranchTraitRateTree(Node* p, std::stringstream& ss);
    void setMeanBranchTraitRates();

    bool isUltrametric();

    // Functions for phenotypic evolution:
    void getPhenotypes(std::string fname);
    void getPhenotypesMissingLatent(std::string fname);

    void  printTraitValues();
    void  initializeTraitValues();
    void  recursiveSetTraitValues(Node* x, double mn, double mx);
    Node* chooseInternalNodeAtRandom();

    void   generateTraitsAllNodesBM(Node* xnode, double varx);
    void   generateTraitsAllNodesFromEventBeta(Node* xnode);
    void   printTraitRange();
    double getTraitMinTip();
    double getTraitMaxTip();

    // speciation-extinction output
    void setMeanBranchSpeciation();
    void setMeanBranchExtinction();

    void echoMeanBranchRates();

    void writeMeanBranchSpeciationTree(Node* p, std::stringstream& ss);
    void writeBranchSpeciationRatesToFile(std::string fname, bool append);
    void writeBranchExtinctionRatesToFile(std::string fname, bool append);


    void writeMeanBranchExtinctionTree(Node* p, std::stringstream& ss);
    void writeNodeSpeciationTree(Node* p, std::stringstream& ss);

    void writeMeanBranchNetDivRateTree(Node* p, std::stringstream& ss);
    void writeBranchPhenotypes(Node* p, std::stringstream& ss);

    // speciation-extinction initialization:

    // File for species-specific values.  
    void initializeSpeciationExtinctionModel(std::string fname);
    void initializeSpeciationExtinctionModel(double sampFraction);
    void printInitialSpeciationExtinctionRates();

    // New fxns for mapping events to nodes:

    double getTotalMapLength();
    void   setTotalMapLength(double x);

    void setCanNodeBeMapped(int ndesc);

    void loadPreviousNodeStates(Tree* ostree);

    // Functions for random access of nodes from temporary nodeset array
    void  setTempInternalNodeArray(Node* p);
    void  tempNodeSetPassDown(Node* p);
    void  clearTempNodeArray();
    Node* getRandomNodeFromTempArray();

    void printNodeBranchRates();

    void computeMeanTraitRatesByNode(Node* x);

    Node* getNodeMRCA(std::string A, std::string B);
    void  passUpFillTempNodeArray(Node* x);
    Node* getNodeByName(std::string A);

    void printNodeTraitRates();

    void printCanHoldEventByNode();

    void echoMeanBranchTraitRates();
};

inline double Tree::getTreeLength()
{
    return treelength;
}


inline Node* Tree::getDownPassNode(int i)
{
    return downPassSeq[i];
}


inline int Tree::getNumberOfNodes()
{
    return (int)downPassSeq.size();
}


inline Node* Tree::getNodeFromDownpassSeq(int i)
{
    return downPassSeq[i];
}


inline void Tree::setStartTime(double x)
{
    _startTime = x;
}


inline void Tree::setTmax(double x)
{
    _tmax = x;
}


inline double Tree::getTmax()
{
    return _tmax;
}


inline Node* Tree::getRoot()
{
    return root;
}


inline void Tree::setExtantStatus(bool x)
{
    _isExtant = x;
}


inline bool Tree::getExtantStatus()
{
    return _isExtant;
}


inline double Tree::getAge()
{
    return _age;
}


inline double Tree::getTotalMapLength()
{
    return _totalMapLength;
}


inline void Tree::setTotalMapLength(double x)
{
    _totalMapLength = x;
}


#endif
