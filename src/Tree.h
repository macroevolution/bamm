#ifndef TREE_H
#define TREE_H

#include "NewickTreeReader.h"
#include "Node.h"

#include <string>
#include <set>
#include <vector>
#include <iosfwd>

class Random;
class Settings;
class branchEvent;
class eventSet;
class Phenotype;
class BranchHistory;
class TraitBranchHistory;
class Node;


// The variance of the root-to-tip lengths must be less than
// this value for the tree to be considered ultrametric
#define ULTRAMETRIC_TOLERANCE 1e-6


// Comparison function used in Tree::setAge()
inline bool compareNodeTime(Node* i, Node* j)
{
    return i->getTime() < j->getTime();
}


class Tree
{

private:

    void setPreOrderNodes(Node* node);
    void setPostOrderNodes(Node* node);

    double calculateTreeLength();
    void setInternalNodes();
    void setNodeTipCounts();

    std::vector<double> terminalPathLengthsToRoot();
    void storeTerminalPathLengthsToRootRecurse
        (Node* node, std::vector<double>& pathLengths);
    bool allValuesAreTheSame(std::vector<double>& list);

    void assertTreeRootBranchLengthIsZero();
    void assertTreeIsBifurcating();
    void assertTreeIsBifurcatingRecurse(Node* node);
    void assertBranchLengthsArePositive();
    void assertBranchLengthsArePositiveRecurse(Node* node);
    void assertTreeIsUltrametric();
    void assertTipsHaveUniqueNames();
    std::vector<std::string> terminalNames();
    void storeTerminalNamesRecurse(Node* node, std::vector<std::string>& names);

    void storeTerminalNodesRecurse(Node* node, std::vector<Node*>& nodes);

    void crossValidateSpecies(const std::vector<std::string>& species);
    void assertIsSubset(const std::vector<std::string>& list1,
                        const std::vector<std::string>& list2,
                        const std::string& list2Name);

    Random& _random;

    Node* root;
    std::vector<Node*> _preOrderNodes;
    std::vector<Node*> _postOrderNodes;

    std::vector<Node*> _internalNodes;

    double _startTime;
    double _tmax;
    bool _isExtant;
    void recursivelyAddNodesToSet(Node* p);
    void rebuildTreeNodeSet();
    double _age; // time to root node, from present

    void setIsLivingTipStatus();
    void setTipStatus();

    // This is the total treelength: also used in mapping events to tree
    double _treeLength;

    // This is a pointer to an object that stores
    // all events that happened.
    eventSet* treeEvents; // What is this doing???

    // New data for mapping events to nodes:
    // Nodes subtended by a branch that can hold an event
    std::set<Node*> mappableNodes;
    double _totalMapLength;

    std::set<Node*> _tempNodeSet;

    NewickTreeReader _treeReader;

public:

    Tree(Random& random, Settings& settings);

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

    int   getNumberOfNodes();
    const std::vector<Node*>& postOrderNodes();

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
    void setRoot(Node* rootNode);

    std::string getNewick();
    void   writeTree(Node* p, std::ostream& ss);

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
    void readTree(const std::string& treeFileName);
    bool isValidChar(char x);

    void setNodeTimes();
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

    void setNodeSpeciationParameters();
    void setNodeExtinctionParameters();
    
    void echoMeanBranchRates();

    void writeMeanBranchSpeciationTree(Node* p, std::stringstream& ss);
    void writeBranchSpeciationRatesToFile(std::string fname, bool append);
    void writeBranchExtinctionRatesToFile(std::string fname, bool append);


    void writeMeanBranchExtinctionTree(Node* p, std::stringstream& ss);
    void writeNodeSpeciationTree(Node* p, std::stringstream& ss);

    void writeMeanBranchNetDivRateTree(Node* p, std::stringstream& ss);
    void writeBranchPhenotypes(Node* p, std::ostream& out);

    // speciation-extinction initialization:

    // File for species-specific values.  
    void initializeSpeciationExtinctionModel(std::string fname);
    void initializeSpeciationExtinctionModel(double sampFraction);
    void printInitialSpeciationExtinctionRates();

    // New fxns for mapping events to nodes:

    double getTotalMapLength();
    void   setTotalMapLength(double x);

    double maxRootToTipLength();

    void setCanNodeBeMapped(int ndesc);

    // Functions for random access of nodes from temporary nodeset array
    void  setTempInternalNodeArray(Node* p);
    void  tempNodeSetPassDown(Node* p);
    void  clearTempNodeArray();
    Node* getRandomNodeFromTempArray();

    Node* getRandomNonRootNode();

    void printNodeBranchRates();

    void computeMeanTraitRatesByNode(Node* x);

    Node* getNodeMRCA(const std::string& A, const std::string& B);
    void  passUpFillTempNodeArray(Node* x);
    Node* getNodeByName(const std::string& A);

    void printNodeTraitRates();

    void echoMeanBranchTraitRates();

    std::vector<Node*> terminalNodes();
    std::vector<double> traitValues();
};


inline double Tree::getTreeLength()
{
    return _treeLength;
}


inline int Tree::getNumberOfNodes()
{
    return (int)_preOrderNodes.size();
}


inline const std::vector<Node*>& Tree::postOrderNodes()
{
    return _postOrderNodes;
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


inline void Tree::setRoot(Node* rootNode)
{
    root = rootNode;
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
