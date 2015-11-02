#include "Random.h"
#include "Settings.h"
#include "Tree.h"
#include "Node.h"
#include "NewickTreeReader.h"
#include "BranchHistory.h"
#include "TraitBranchEvent.h"
#include "Log.h"
#include "Stat.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>


Tree::Tree(Random& random, Settings& settings) : _random(random)
{
    readTree(settings.get("treefile"));

    setPreOrderNodes(root);
    setPostOrderNodes(root);

    setStartTime(0.0);
    setNodeTimes();
    setAge();

    // Check tree integrity
    assertTreeRootBranchLengthIsZero();
    assertTreeIsBifurcating();
    assertBranchLengthsArePositive();
    if (settings.get<bool>("checkUltrametric")) {
        assertTreeIsUltrametric();
    }
    assertTipsHaveUniqueNames();

    // Output stuff here
    log() << "Tree contains " << getNumberTips() << " taxa.\n";

    _totalMapLength = 0.0;

    setBranchingTimes(root);

    _treeLength = calculateTreeLength();

    setInternalNodes();
    setNodeTipCounts();

    // Initialize tree according to model type
    // TODO: This should be handled in a better way
    if (settings.get("modeltype") == "speciationextinction") {
        if (settings.get<bool>("useGlobalSamplingProbability")) {
            initializeSpeciationExtinctionModel
                (settings.get<double>("globalSamplingFraction"));
        } else {
            initializeSpeciationExtinctionModel
                (settings.get("sampleProbsFilename"));
        }

        setCanNodeHoldEventByDescCount
            (settings.get<int>("minCladeSizeForShift"));
        setTreeMap(getRoot());
    } else if (settings.get("modeltype") == "trait") {
        setAllNodesCanHoldEvent();
        setTreeMap(getRoot());
        getPhenotypesMissingLatent(settings.get("traitfile"));
        initializeTraitValues();
    }
}


double Tree::calculateTreeLength()
{
    double treeLength = 0.0;

    for (int i = 0; i < (int)_preOrderNodes.size(); ++i) {
        treeLength += _preOrderNodes[i]->getBrlen();
    }

    return treeLength;
}


void Tree::setInternalNodes()
{
    for (int i = 0; i < (int)_preOrderNodes.size(); ++i) {
        if (_preOrderNodes[i]->isInternal()) {
            _internalNodes.push_back(_preOrderNodes[i]);
        }
    }
}


void Tree::setNodeTipCounts()
{
    for (int i = 0; i < (int)_preOrderNodes.size(); ++i) {
        int count = getDescTipCount(_preOrderNodes[i]);
        _preOrderNodes[i]->setTipDescCount(count);
    }
}


void Tree::setPreOrderNodes(Node* node)
{
    if (node != NULL) {
        _preOrderNodes.push_back(node);
        setPreOrderNodes(node->getLfDesc());
        setPreOrderNodes(node->getRtDesc());
    }
}


void Tree::setPostOrderNodes(Node* node)
{
    if (node != NULL) {
        setPostOrderNodes(node->getLfDesc());
        setPostOrderNodes(node->getRtDesc());
        _postOrderNodes.push_back(node);
    }
}


Tree::~Tree()
{
    for (int i = 0; i < (int)_preOrderNodes.size(); ++i) {
        delete _preOrderNodes[i];
    }
}

/* Tree map goes from low to high, tipwise (e.g.,
 smaller values closer to root of tree
 high values at TIP of tree

 MapStart is at the base of the branch (to the root)
 MapEnd is at node itself.

  // Deprecated
 */
/*
void Tree::setTreeMap(Node* p){

    p->setMapStart(treelength);
    treelength += p->getBrlen();
    p->setMapEnd(treelength);

    if (p->getRtDesc() != NULL){
        setTreeMap(p->getRtDesc());
    }
    if (p->getLfDesc() != NULL){
        setTreeMap(p->getLfDesc());
    }

}
*/

/*
 setTreeMap
 Requires bool _canHoldEvent attribute on nodes
 allows you to exclude nodes that fail to meet some criterion,
 e.g., no "events" on terminal branches.

 */

void Tree::setTreeMap(Node* p)
{
    p->setMapStart(_totalMapLength);
    _totalMapLength += p->getBrlen();
    p->setMapEnd(_totalMapLength);

    if (p->getRtDesc() != NULL) {
        if (p->getRtDesc()->getCanHoldEvent()) {
            setTreeMap(p->getRtDesc());
        }
    }
    if (p->getLfDesc() != NULL) {
        if (p->getLfDesc()->getCanHoldEvent()) {
            setTreeMap(p->getLfDesc());
        }
    }
}

/*

 Function to recover absolute time (0 at root, T at present)
 from map time. Critical for time-homogeneous birth-death model

 */

// Should NEVER be applied to value of 0.0 (eg at the root).
double Tree::getAbsoluteTimeFromMapTime(double x)
{
    double abstime = 0.0;
    bool done = false;
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if (x >= (*i)->getMapStart()  && x < (*i)->getMapEnd()) {
            double delta = x - (*i)->getMapStart(); // difference in times...
            abstime = (*i)->getTime() - delta;
            done = true;
        }
    }
    if (done == false) {
        std::cout << "could not find abs time from map time \n";
        std::cout << "Tree::getAbsoluteTimeFromMapTime() " << std::endl;
        throw;
    }
    return abstime;
}


// Get number of descendant nodes from a given node
int Tree::getDescNodeCount(Node* p)
{
    int count = 0;
    if (p->getLfDesc() != NULL) {
        count++;
        count += getDescNodeCount(p->getLfDesc());
    }
    if (p->getRtDesc() != NULL) {
        count++;
        count += getDescNodeCount(p->getRtDesc());
    }
    return count;
}

/*
 mapEventToTree
 Event is mapped to tree by map value; each "mappable" branch
 on tree has start and end values for mapping that define a unique interval
 of a real number line.
 */

Node* Tree::mapEventToTree(double x)
{
    Node* y = NULL;
    for (std::set<Node*>::iterator i = mappableNodes.begin();
            i != mappableNodes.end(); i++) {
        if ( x > (*i)->getMapStart() && x <= (*i)->getMapEnd()) {
            y = (*i);
        }
    }
    if (y == NULL) {
        std::cout << "error: unmapped event\n" << std::endl;
        std::cout << "position: " << x << std::endl;
    }
    //std::cout << "y in map: " << y << std::endl;
    return y;
}


void Tree::printNodeMap()
{
    for (std::set<Node*>::iterator i = mappableNodes.begin();
            i != mappableNodes.end(); i++) {
        std::cout << (*i) << "\t" << (*i)->getAnc() << "\t" << (*i)->getMapStart() << "\t"
             << (*i)->getMapEnd() << std::endl;
    }
}


// This requires that p be an internal node to begin with!
void Tree::setTempInternalNodeArray(Node* p)
{
    if (p->getRtDesc() == NULL && p->getLfDesc() == NULL) {
        std::cout << "Problem: sent terminal node to setTempInternalNodeArray" << std::endl;
        throw;
    }
    tempNodeSetPassDown(p);
}


void Tree::tempNodeSetPassDown(Node* p)
{
    if (p->getLfDesc() != NULL && p->getRtDesc() != NULL) {
        _tempNodeSet.insert(p);
        tempNodeSetPassDown(p->getLfDesc());
        tempNodeSetPassDown(p->getRtDesc());
    }
}


void Tree::clearTempNodeArray()
{
    _tempNodeSet.clear();
}


Node* Tree::getRandomNodeFromTempArray()
{
    int chosen = _random.uniformInteger(0, (int)_tempNodeSet.size() - 1);
    int myit = 0;
    Node* xnode = (*_tempNodeSet.begin());

    for (std::set<Node*>::iterator i = _tempNodeSet.begin();
            i != _tempNodeSet.end(); i++) {
        if (myit == chosen) {
            xnode = (*i);
        }
        myit++;
    }
    return xnode;
}


Node* Tree::getRandomNonRootNode()
{
    // Start at index = 1 because the root is at index = 0
    int randomIndex = _random.uniformInteger(1, _preOrderNodes.size() - 1);
    return _preOrderNodes[randomIndex];
}


void Tree::setBranchLengths()
{
    // Set brlens:
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i) == root) {
            double timeX = (*i)->getTime() - _startTime;
            (*i)->setBrlen( timeX );
        } else {
            double timeX = (*i)->getTime() - ( (*i)->getAnc() )->getTime();
            (*i)->setBrlen( timeX );
        }
    }
}


std::string Tree::getNewick()
{
    std::stringstream ss;
    writeTree(root, ss);
    std::string newick = ss.str();
    newick.append(";");
    return newick;
}


void Tree::writeTree(Node* p, std::ostream& ss)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "") {
            ss << p->getIndex() << ":" << p->getBrlen();
        } else {
            ss << p->getName() << ":" << p->getBrlen();
        }
    } else {
        ss << "(";
        writeTree(p->getLfDesc(), ss);
        ss << ",";
        writeTree(p->getRtDesc(), ss);
        ss << "):" << p->getBrlen();
    }
}


int Tree::getNumberTips()
{
    int count = 0;
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ( (*i)->getDescCount() == 0 ) {
            count++;
        }
    }
    return count;
}


void Tree::setTipStatus()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if (((*i)->getLfDesc() == NULL) && ((*i)->getRtDesc() == NULL )) {
            (*i)->setIsTip(true);
        } else {
            (*i)->setIsTip(false);
        }
    }
}


int Tree::getNumberExtantTips()
{
    setExtantStatus();
    int ntaxa = 0;
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if (((*i)->getLfDesc() == NULL) && ((*i)->getExtantStatus() == 1) ) {
            ntaxa++;
        }
    }
    return ntaxa;
}



/* setExtantStatus
        value = 1 if node is a leaf AND if
        node is extant
        OR if node is internal

        extinct leaves get value = 0

*/


void Tree::setExtantStatus()
{
    double tx = root->getTime();
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getTime() > tx) {
            tx = (*i)->getTime();
        }
    }
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getLfDesc() != NULL && (*i)->getRtDesc() != NULL) {
            (*i)->setExtantStatus(1);
        } else {
            if ((*i)->getTime() > (tx - 0.000001) ) {
                (*i)->setExtantStatus(1);
            } else {
                (*i)->setExtantStatus(0);
            }
        }
        //std::cout << "Node\t" << (*i)->getIndex() << "\t" << (*i)->getExtantStatus() << std::endl;
    }
}


/*
 This function iterates over all nodes.
 If the node is a tip, and ALIVE, it gets flagged with 1. Otherwise 0.
 This can be distinguished from setExtantStatus,
    which flags status of internal nodes as well.
 */
void Tree::setIsLivingTipStatus()
{
    double TOL = 0.000001;

    // Get maximum age.
    double tx = root->getTime();
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getTime() > tx) {
            tx = (*i)->getTime();
        }
    }

    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL) {
            double test = abs((*i)->getTime() - tx );
            if (test < TOL) {
                (*i)->setIsLivingTip(true);
            } else {
                (*i)->setIsLivingTip(false);
            }
        } else {
            (*i)->setIsLivingTip(false);
        }
    }
}


void Tree::rebuildTreeNodeSet()
{
    _preOrderNodes.clear();
    recursivelyAddNodesToSet(root);
}


void Tree::recursivelyAddNodesToSet(Node* p)
{
    if (p != NULL) {
        _preOrderNodes.push_back(p);
        recursivelyAddNodesToSet(p->getLfDesc());
        recursivelyAddNodesToSet(p->getRtDesc());
    }
}



/* Function to drop tips from a tree
  where ALL TIPS are extant.

 You must call pruneExtinctTaxa() on the tree
  before using this function, or results will be
  garbage.

*/


void Tree::pruneExtinctTaxa()
{
    //assume most recent time is extant
    setExtantStatus();
    fixExtinct(root);

    rebuildTreeNodeSet();

    setBranchLengths();
}


/*

void Tree::fixExtinct(Node * p)

A recursive function for eliminating extinct taxa
 and all 'cherry' nodes created by eliminating them

 *should* generate ultrametric tree

*/


void Tree::fixExtinct(Node* p)
{
    if (p->getDescCount() == 0) {
        std::cout << "terminal node: should never get here\n" << std::endl;
    }
    // separate recursion from modification:
    // Recursion step:
    if (p->getLfDesc()->getDescCount() == 2) {
        fixExtinct(p->getLfDesc());
    }
    if (p->getRtDesc()->getDescCount() == 2) {
        fixExtinct(p->getRtDesc());
    }

    // modification step:

    if (p->getLfDesc()->getExtantStatus() == 0 &&
            p->getRtDesc()->getExtantStatus() == 1) {
        //std::cout << "here" << std::endl;
        if (p == root) {
            root = p->getRtDesc();
            p->getRtDesc()->nullifyAnc();
            p->nullifyLfDesc();
            p->nullifyRtDesc();
            //nodes.erase(p);
            delete p;
        } else if (p->getAnc()->getLfDesc() == p) {
            p->getAnc()->setLfDesc(p->getRtDesc());
            p->getRtDesc()->setAnc(p->getAnc());
            //nodes.erase(p);
            delete p;
        } else if (p->getAnc()->getRtDesc() == p) {
            p->getAnc()->setRtDesc(p->getRtDesc());
            p->getRtDesc()->setAnc(p->getAnc());
            //nodes.erase(p);
            delete p;
        } else {
            std::cout << "problem: p state invalid" << std::endl;
        }

    } else if (p->getLfDesc()->getExtantStatus() == 1 &&
               p->getRtDesc()->getExtantStatus() == 0) {

        if (p == root) {
            //std::cout << "p is root" << std::endl;
            //std::cout << root << "\t" << p << "\t" << std::endl;
            root = p->getLfDesc();
            //std::cout << p->getLfDesc() << "\t" << root << std::endl;
            p->getLfDesc()->nullifyAnc();
            p->nullifyLfDesc();
            p->nullifyRtDesc();
            //nodes.erase(p);
            delete p;
        } else if (p->getAnc()->getLfDesc() == p) {
            p->getAnc()->setLfDesc(p->getLfDesc());
            p->getLfDesc()->setAnc(p->getAnc());
            //nodes.erase(p);
            delete p;
        } else if (p->getAnc()->getRtDesc() == p) {
            p->getAnc()->setRtDesc(p->getLfDesc());
            p->getLfDesc()->setAnc(p->getAnc());
            //nodes.erase(p);
            delete p;
        } else {
            std::cout << "invalid p state" << std::endl;
        }

    } else if (p->getLfDesc()->getExtantStatus() == 0 &&
               p->getRtDesc()->getExtantStatus() == 0) {
        if (p == root) {
            std::cout << "Problem: tree extinct" << std::endl;
        } else {
            // collapse nodes
            //nodes.erase(p->getRtDesc());
            //nodes.erase(p->getLfDesc());
            delete p->getRtDesc();
            delete p->getLfDesc();
            p->nullifyLfDesc();
            p->nullifyRtDesc();
            p->setExtantStatus(0);
            //fixExtinct(p->getAnc());
        }

    } else if (p->getLfDesc()->getExtantStatus() == 1 &&
               p->getRtDesc()->getExtantStatus() == 1) {
        //continue
    } else {
        std::cout << "problem: invalid condition" << std::endl;
    }
}


/* deleteExtinctNodes
 This function does:
    1. find extinct node
    2. set pointer pointing to node to NULL
    3. delete node

*/


void Tree::deleteExtinctNodes()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        // if node is leaf:
        std::cout << (*i)->getIndex() << std::endl;
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL) {
            // if nodes is extinct
            if ((*i)->getExtantStatus() == 0) {
                if ((*i)->getAnc()->getLfDesc() == (*i)) {
                    (*i)->getAnc()->setLfDesc(NULL);
                    //delete (*i);
                } else if ((*i)->getAnc()->getRtDesc() == (*i) ) {
                    (*i)->getAnc()->setRtDesc(NULL);
                    //delete (*i);
                } else {
                    std::cout << "problem in deleteExtinctNodes()" << std::endl;
                }
            }
        }
    }
}


void Tree::readTree(const std::string& treeFileName)
{
    std::ifstream treeFileStream(treeFileName.c_str());

    log() << "\nReading tree from file <" << treeFileName << ">.\n";

    if (!treeFileStream.good()) {
        log(Error) << "Invalid file name for phylogenetic tree\n";
        std::exit(1);
    }

    _treeReader.read(treeFileStream, *this);
}


void Tree::setNodeTimes()
{
    // Handle the root node (the first element in a pre-order traversal)
    if (getNumberOfNodes() > 0) {
        _preOrderNodes[0]->setTime(0.0);
    }

    for (int i = 1; i < (int)_preOrderNodes.size(); ++i) {
        Node* node = _preOrderNodes[i];
        node->setTime(node->getBrlen() + node->getAnc()->getTime());
    }
}


bool Tree::isValidChar(char x)
{
    if (x == ')' || (x == '(') || (x == ':') || (x == ',') || (x == ';')) {
        return false;
    } else {
        return true;
    }
}


void Tree::setAge()
{
    Node* nodeWithMaxTime = *std::max_element
        (_preOrderNodes.begin(), _preOrderNodes.end(), compareNodeTime);
    _age = nodeWithMaxTime->getTime();
}


void Tree::setBranchingTimes(Node* p)
{
    double bt = getAge() - p->getTime();
    p->setBranchTime(bt);
    if (p->getLfDesc() != NULL && p->getRtDesc() != NULL) {
        setBranchingTimes(p->getLfDesc());
        setBranchingTimes(p->getRtDesc());
    }

}


std::vector<double>  Tree::getBranchingTimes()
{
    double TOL = 0.0000001;
    std::vector<double> btimes;

    //std::cout << "age :" << _age << std::endl;

    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); i++) {
        double tmp = _age - (*i)->getTime();
        //std::cout << "Time: " << (*i)->getTime() << "\t" << tmp << std::endl;

        if (tmp > TOL) {
            btimes.push_back(tmp);
        }
    }
    std::sort(btimes.begin(), btimes.end() );
    std::reverse(btimes.begin(), btimes.end());
    return btimes;
}


void Tree::writeMeanBranchTraitRateTree(Node* p, std::stringstream& ss)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "") {
            ss << p->getIndex() << ":" << p->getMeanBeta();
        } else {
            ss << p->getName() << ":" << p->getMeanBeta();
        }
    } else {
        ss << "(";
        writeMeanBranchTraitRateTree(p->getLfDesc(), ss);
        ss << ",";
        writeMeanBranchTraitRateTree(p->getRtDesc(), ss);
        ss << "):" << p->getMeanBeta();
    }
}


/*

 Set mean speciation rates on each branch by going over all nodes
 and accessing BranchHistory attribute.

 If multiple events on a particular branch, the meanBranchRate is
 just the arithmetic average of the rates (for now, anyway).

 */

void Tree::setNodeSpeciationParameters()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i){
        (*i)->computeAndSetNodeSpeciationParams();
    }
}

void Tree::setNodeExtinctionParameters()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i){
        (*i)->computeAndSetNodeExtinctionParams();
    }
}


// Update both mean speciation rates on branch in addition to node speciation rate.
void Tree::setMeanBranchSpeciation()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        (*i)->computeNodeBranchSpeciationParams();
    }
}

// this will update extinction for both branch mean rates as well as
// node-associated rates.

void Tree::setMeanBranchExtinction()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        (*i)->computeNodeBranchExtinctionParams();
    }
}


void Tree::echoMeanBranchRates()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        Node* x = (*i);
        std::cout << x << "\t" << x->getMeanSpeciationRate() << "\t" <<
             x->getMeanExtinctionRate() << "\t" << x->getNodeLambda() << std::endl;
    }
}


void Tree::echoMeanBranchTraitRates()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        Node* x = (*i);
        std::cout << x->getName() << "\t" << x->getMeanBeta() << std::endl;
    }
}


void Tree::writeBranchSpeciationRatesToFile(std::string fname, bool append)
{
    std::stringstream outdata;

    writeMeanBranchSpeciationTree(root, outdata);

    outdata << ";";

    std::ofstream outStream;
    if (append == true) {
        outStream.open(fname.c_str(), std::ofstream::app);
    } else {
        outStream.open(fname.c_str(), std::ofstream::trunc);
    }
    outStream << outdata.str() << std::endl;
    outStream.close();
}


void Tree::writeBranchExtinctionRatesToFile(std::string fname, bool append)
{
    std::stringstream outdata;
    writeMeanBranchExtinctionTree(root, outdata);
    outdata << ";";

    std::ofstream outStream;
    if (append == true) {
        outStream.open(fname.c_str(), std::ofstream::app);
    } else {
        outStream.open(fname.c_str(), std::ofstream::trunc);
    }
    outStream << outdata.str() << std::endl;
    outStream.close();
}


void Tree::writeMeanBranchSpeciationTree(Node* p, std::stringstream& ss)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "") {
            ss << p->getIndex() << ":" << p->getMeanSpeciationRate();
        } else {
            ss << p->getName() << ":" << p->getMeanSpeciationRate();
        }
    } else {
        ss << "(";
        writeMeanBranchSpeciationTree(p->getLfDesc(), ss);
        ss << ",";
        writeMeanBranchSpeciationTree(p->getRtDesc(), ss);
        ss << "):" << p->getMeanSpeciationRate();
    }
}


void Tree::writeNodeSpeciationTree(Node* p, std::stringstream& ss)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "") {
            ss << p->getIndex() << ":" << p->getNodeLambda();
        } else {
            ss << p->getName() << ":" << p->getNodeLambda();
        }
    } else {
        ss << "(";
        writeNodeSpeciationTree(p->getLfDesc(), ss);
        ss << ",";
        writeNodeSpeciationTree(p->getRtDesc(), ss);
        ss << "):" << p->getNodeLambda();
    }
}


void Tree::writeMeanBranchExtinctionTree(Node* p, std::stringstream& ss)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "") {
            ss << p->getIndex() << ":" << p->getMeanExtinctionRate();
        } else {
            ss << p->getName() << ":" << p->getMeanExtinctionRate();
        }
    } else {
        ss << "(";
        writeMeanBranchExtinctionTree(p->getLfDesc(), ss);
        ss << ",";
        writeMeanBranchExtinctionTree(p->getRtDesc(), ss);
        ss << "):" << p->getMeanExtinctionRate();
    }
}


void Tree::writeMeanBranchNetDivRateTree(Node* p, std::stringstream& ss)
{

    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "") {
            ss << p->getIndex() << ":" << (p->getMeanSpeciationRate() -
                                           p->getMeanExtinctionRate());
        } else {
            ss << p->getName() << ":" << (p->getMeanSpeciationRate() -
                                          p->getMeanExtinctionRate());
        }
    } else {
        ss << "(";
        writeMeanBranchNetDivRateTree(p->getLfDesc(), ss);
        ss << ",";
        writeMeanBranchNetDivRateTree(p->getRtDesc(), ss);
        ss << "):" << (p->getMeanSpeciationRate() - p->getMeanExtinctionRate());
    }
}


void Tree::writeBranchPhenotypes(Node* p, std::ostream& out)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "") {
            out << p->getIndex() << ":" << p->getTraitValue();
        } else {
            out << p->getName() << ":" << p->getTraitValue();
        }
    } else {
        out << "(";
        writeBranchPhenotypes(p->getLfDesc(), out);
        out << ",";
        writeBranchPhenotypes(p->getRtDesc(), out);
        out << "):" << p->getTraitValue();
    }
}


/*
Tree::getPhenotypes
Read file. First column = species name exactly as matching in phylogeny.
 Second column, tab-delimited: phenotype value on appropriate scale
    (e.g., already log-transformed).

 */

void Tree::getPhenotypes(std::string fname)
{
    std::ifstream infile(fname.c_str());
    log() << "\nReading phenotypes from file <" << fname.c_str() << ">\n";
    std::vector<std::string> stringvec;
    std::vector<std::string> spnames;
    std::vector<double> traits;

    if (!infile.good()) {
        log(Error) << "Error reading file.\n";
        std::exit(1);
    }

    while (infile) {
        std::string tempstring;
        getline(infile, tempstring, '\t');
        //std::cout << tempstring << "\n" << std::endl;
        spnames.push_back(tempstring);
        getline(infile, tempstring, '\n');
        traits.push_back(atof(tempstring.c_str()));

        // this OK?
        if (infile.peek() == EOF) {
            break;
        }
    }

    infile.close();

    log() << "Read " << traits.size() << " species with trait data\n";

    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL ) {
            for (std::vector<std::string>::size_type k = 0; k < spnames.size(); k++) {
                if ((*i)->getName() == spnames[k]) {
                    (*i)->setTraitValue(traits[k]);
                    (*i)->setIsTraitFixed(true);
                }
            }
            if ((*i)->getIsTraitFixed() == false) {
                std::cout << "error - failed to set a terminal state\n" << std::endl;
            }
        } else {
            (*i)->setTraitValue(0);
            (*i)->setIsTraitFixed(false);
        }
    }
}


// This and the function above could be combined into one -- JWB
void Tree::getPhenotypesMissingLatent(std::string fileName)
{
    std::ifstream inputFile(fileName.c_str());

    if (!inputFile) {
        log(Error) << "Could not read trait values from file "
            << "<<" << fileName << ">>.\n";
        std::exit(1);
    }

    log() << "Reading traits from file <<" << fileName << ">>.\n";

    std::vector<std::string> speciesNames;
    std::vector<double> traitValues;

    std::string speciesName;
    double traitValue;

    while (inputFile >> speciesName) {
        if (!(inputFile >> traitValue)) {
            exitWithError("Invalid trait value for <" + speciesName + ">.");
        }

        speciesNames.push_back(speciesName);
        traitValues.push_back(traitValue);
    }

    inputFile.close();

    log() << "Read " << traitValues.size() << " species with trait data.\n";

    int missingTerminalCount = 0;

    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL ) {
            for (int k = 0; k < (int)speciesNames.size(); k++) {
                if ((*i)->getName() == speciesNames[k]) {
                    (*i)->setTraitValue(traitValues[k]);
                    (*i)->setIsTraitFixed(true);
                }
            }
            if ((*i)->getIsTraitFixed() == false) {
                missingTerminalCount++;
            }
        } else {
            (*i)->setTraitValue(0);
            (*i)->setIsTraitFixed(false);
        }
    }

    if (missingTerminalCount > 0) {
        log(Warning) << "Missing data for < " << missingTerminalCount
                     << " > species.\n"
                     << "These will be treated as latent variables in "
                     << "analysis\n";
    }
}


// recursive function to generate trait values...
void Tree::generateTraitsAllNodesBM(Node* xnode, double varx)
{
    if (xnode->getLfDesc() != NULL && xnode->getRtDesc() != NULL) {
        // do left:
        double vx = xnode->getLfDesc()->getBrlen() * varx;
        double newTrait = xnode->getTraitValue() +
            _random.normal(0.0, std::sqrt(vx));
        xnode->getLfDesc()->setTraitValue(newTrait);
        generateTraitsAllNodesBM(xnode->getLfDesc(), varx);

        vx = xnode->getRtDesc()->getBrlen() * varx;
        newTrait = xnode->getTraitValue() +
            _random.normal(0.0, std::sqrt(vx));
        xnode->getRtDesc()->setTraitValue(newTrait);
        generateTraitsAllNodesBM(xnode->getRtDesc(), varx);


    } else if (xnode->getLfDesc() == NULL && xnode->getRtDesc() == NULL) {
        xnode->setIsTraitFixed(true);
    } else {
        std::cout << "error in Tree::generateTraitsAllNodesBM" << std::endl;
        throw;
    }
    //std::cout << xnode->getTraitValue() << std::endl;
}


/*
    chooseInternalNodeAtRandom()
        have checked this to make sure distribution of sampled nodes is uniform
        appears to be fine.

 */

Node* Tree::chooseInternalNodeAtRandom()
{
    int snode = _random.uniformInteger(0, (int)_internalNodes.size() - 1);
    std::vector<Node*>::iterator myIt = _internalNodes.begin();

    for (int i = 0; i < snode; i++ ) {
        myIt++;
    }
    return  (*myIt);
}



/*
    This sets all internal nodes equal to the estimated sampling fraction for that particular node.
    Ei values for all tip nodes get set;
    Di values for all tip nodes get set.
    Etip gets set for all nodes for global sampling probability

*/

// Modified 11.23.2014 to set D0 and E0 probabilities
//    for fossil tips (probable extinct lineages but tips)

// TODO: check initialization for incomplete sampling for fossil process
void Tree::initializeSpeciationExtinctionModel(double sampFrac)
{
    double  speciationInit = sampFrac;
    double  extinctionInit = (double)1 - sampFrac;

    
    
    for (std::vector<Node*>::iterator myIt = _preOrderNodes.begin();
            myIt != _preOrderNodes.end(); ++myIt) {
        if ((*myIt)->getLfDesc() == NULL && (*myIt)->getRtDesc() == NULL) {

            (*myIt)->setDinit(speciationInit);
            (*myIt)->setEinit(extinctionInit);
            
            // This line unnecessary
            //bool isExtant = (std::abs(getAge() - (*myIt)->getTime() )) <= 0.0001;
          
            // FOSSIL: initial conditions should be same as for
            //     extant tip
            
            // These lines are bad logic & bad math I think
            //if (isExtant){
            //  (*myIt)->setDinit(speciationInit);
            //  (*myIt)->setEinit(extinctionInit);
            //}else{
            //  (*myIt)->setDinit(extinctionInit);
            //  (*myIt)->setEinit(speciationInit);
            //}
            

        }
        (*myIt)->setEtip(extinctionInit); // Set
    }
}

// This does not work with the fossil process yet.

void Tree::initializeSpeciationExtinctionModel(std::string fname)
{
    
    // Tree must be ultrametric to use this option
    // TODO: this check for ultrametric must be more informative
    
    assertTreeIsUltrametric();
    
    std::ifstream infile(fname.c_str());

    
    
    if (!infile.good()) {
        log(Error) << "Bad sampling fraction file.\n";
        std::exit(1);
    }

    log() << "Reading sampling fractions from file <<" << fname << ">>...\n";

    std::vector<std::string> stringvec;
    std::vector<std::string> spnames;
    std::vector<std::string> spfamilies;
    std::vector<double> sfracs;

    std::string tempstring;

    // First number in file is sampling probability for "backbone" of the tree
    getline(infile, tempstring, '\n');
    double backboneSampProb = std::atof(tempstring.c_str());
    double backboneInitial = 1.0 - backboneSampProb;

    // std::atof returns 0.0 if the conversion fails
    if (backboneSampProb == 0.0) {
        log(Error) << "The first line of the sampling probability file\n"
                   << "must be the global sampling probability (> 0.0).\n";
        std::exit(1);
    }

    while (infile) {
        getline(infile, tempstring);
        std::istringstream inStream(tempstring);

        std::string spname;
        std::string spfamily;
        double sfrac = 0.0;

        bool eof = false;

        eof = inStream.eof();
        inStream >> spname;

        eof = inStream.eof();
        inStream >> spfamily;

        eof = inStream.eof();
        inStream >> sfrac;

        if (eof) {
            log(Error) << "Sampling probability file is not formatted "
                       << "properly.\nPlease see the documentation.\n";
            std::exit(1);
        }

        if (sfrac <= 0.0 || sfrac > 1.0) {
            log(Error) << "In line with species <<" << spname << ">> and "
                       << "family name\n<<" << spfamily << ">>, "
                       << "sampling fraction must be greater than 0 and less\n"
                       << "than or equal to 1. This error may also occur if\n"
                       << "the input line is not formatted properly.\n";
            std::exit(1);
        }

        spnames.push_back(spname);
        spfamilies.push_back(spfamily);
        sfracs.push_back(sfrac);

        if (infile.peek() == EOF) {
            break;
        }
    }

    infile.close();

    std::cout << "Read a total of " << sfracs.size() << " initial values.\n";

    crossValidateSpecies(spnames);

    int counter = 0;
    for (std::vector<Node*>::iterator i = _postOrderNodes.begin();
            i != _postOrderNodes.end(); i++) {

        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL ) {
            for (std::vector<std::string>::size_type k = 0; k < spnames.size(); k++) {
                if ((*i)->getName() == spnames[k]) {
                    double Einit = (double)1 - sfracs[k];
                    double Dinit = sfracs[k];

                    (*i)->setEinit(Einit);
                    (*i)->setEtip(Einit);
                    (*i)->setDinit(Dinit);
                    (*i)->setCladeName(spfamilies[k]);
                    counter++;
                }
                //std::cout << spfamilies[k] << std::endl;
            }

            if ((*i)->getEinit() == -1) {
                log(Warning) << "The species " << (*i)->getName() << " "
                    << "has an E_init value of -1.\n"
                    << "Check that this species is spelled correctly "
                    << "in the sampling file.\n";
            }
        } else {
            // Node is internal
            if ((*i)->getLfDesc()->getCladeName() == (*i)->getRtDesc()->getCladeName()) {
                // node *i belongs to same clade and inherits their sampling probability:
                double sprob = (*i)->getLfDesc()->getEtip();
                (*i)->setEtip(sprob);
                (*i)->setCladeName((*i)->getLfDesc()->getCladeName());
            } else {
                std::string cname = "backbone";
                (*i)->setCladeName(cname);
                (*i)->setEtip(backboneInitial);
            }
        }

        if ((*i)->getCladeName() == "") {
            log(Error) << "There are unset clade names.\n";
            std::exit(1);
        }
    }

    int tcount = 0;
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        Node* x = (*i);
        if (x->getEtip() < 0) {
            tcount++;
        }
    }
    std::cout << "Set a total of < " << counter << " > tips nodes for Ei & Di" << std::endl;
    std::cout << "Failed to set < " << tcount << " > internal node eTip values" << std::endl;
}


// Assert that all elements in species are in the tree and vice versa
void Tree::crossValidateSpecies(const std::vector<std::string>& species)
{
    const std::vector<std::string>& treeSpecies = terminalNames();

    assertIsSubset(species, treeSpecies, "tree");
    assertIsSubset(treeSpecies, species, "file");
}


void Tree::assertIsSubset(const std::vector<std::string>& list1,
                          const std::vector<std::string>& list2,
                          const std::string& list2Name)
{
    std::vector<std::string>::const_iterator it;
    for (it = list1.begin(); it != list1.end(); ++it) {
        if (std::find(list2.begin(), list2.end(), *it) == list2.end()) {
            log(Error) << "<<" << *it << ">> is not in " << list2Name << ".\n";
            std::exit(1);
        }
    }
}


void Tree::printInitialSpeciationExtinctionRates()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        std::cout << (*i) << "\t" << (*i)->getDinit() << "\t"
            << (*i)->getEinit()  << "\t" << (*i)->getBrlen() << std::endl;
    }
}


/*

 setCanNodeBeMapped

 A node can be mapped IFF:
    1. It contains a min of ndesc tip descendants (including itself)
        with valid tip data (e.g., node->getIsNodeFixed == true)

 */

void Tree::setCanNodeBeMapped(int ndesc)
{
    for (std::vector<Node*>::iterator myIt = _preOrderNodes.begin();
            myIt != _preOrderNodes.end(); ++myIt) {
        int dcount = countDescendantsWithValidTraitData((*myIt));
        if (dcount >= ndesc) {
            (*myIt)->setCanHoldEvent(true);
            mappableNodes.insert((*myIt));
        }
    }
    std::cout << "Number of mappable nodes: < " << mappableNodes.size() << " >" << std::endl;
}


int Tree::getDescTipCount(Node* p)
{
    int count = 0;
    if ((p->getLfDesc() == NULL) & (p->getRtDesc() == NULL)) {
        count++;
    } else {
        count += getDescTipCount(p->getLfDesc());
        count += getDescTipCount(p->getRtDesc());
    }
    return count;
}

/*
    counts number of descendant tips from a given node,
    subject to the condition that all tips have valid
    trait data, e.g., no missing values...

 */
int Tree::countDescendantsWithValidTraitData(Node* p)
{
    int count = 0;
    if ((p->getLfDesc() == NULL) && (p->getRtDesc() == NULL)) {
        if (p->getIsTraitFixed())
            count++;
    } else {
        count += countDescendantsWithValidTraitData(p->getLfDesc());
        count += countDescendantsWithValidTraitData(p->getRtDesc());
    }
    return count;
}

/*

 Sets status for all nodes such that
 _canHoldEvent = true
 This is only used for the 'speciation/extinction' version of BAMM
 where no trait data are necessary.

 */

void Tree::setAllNodesCanHoldEvent()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        (*i)->setCanHoldEvent(true);
        mappableNodes.insert((*i));
    }
}


void Tree::setCanNodeHoldEventByDescCount(int x)
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getTipDescCount() >= x) {
            (*i)->setCanHoldEvent(true);
            mappableNodes.insert((*i));
        }
    }
}


void Tree::setTreeMap(int ndesc)
{
    // set bool flag on each node, depending on whether event can be mapped...
    setCanNodeBeMapped(ndesc);
    setTreeMap(root);
    std::cout << "Map length: " << getTotalMapLength() << std::endl;
    std::cout << "Total mappable nodes: " << mappableNodes.size() << std::endl;
}


void Tree::printTraitRange()
{
    double mx = root->getTraitValue();
    double mn = root->getTraitValue();
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        double ctrait = (*i)->getTraitValue();
        if (ctrait < mn) {
            mn = ctrait;
        }
        if (ctrait > mx) {
            mx = ctrait;
        }
    }
    std::cout << "Min trait value < " << mn << " >\tMax value < " << mx << " >" << std::endl;
}


double Tree::getTraitMaxTip()
{
    double ctrait = root->getTraitValue();
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getLfDesc() == NULL && (*i)->getTraitValue() > ctrait ) {
            ctrait = (*i)->getTraitValue();
        }
    }
    return ctrait;
}


double Tree::getTraitMinTip()
{
    double ctrait = root->getTraitValue();
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getLfDesc() == NULL && (*i)->getTraitValue() < ctrait )
            ctrait = (*i)->getTraitValue();
    }
    return ctrait;
}


void Tree::printNodeBranchRates()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        std::cout << (*i)->getMeanSpeciationRate() << "\t" << (*i)->getNodeLambda() << std::endl;
    }
}


void Tree::printNodeTraitRates()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        std::cout << (*i) << "\t" << (*i)->getMeanBeta() << std::endl;
    }
}


void Tree::setMeanBranchTraitRates()
{
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        computeMeanTraitRatesByNode((*i));
    }
}


/*
    This function is replicated for speciation and extinction as part of
    class node - it seems more efficient to put it here with class tree.


*/

void Tree::computeMeanTraitRatesByNode(Node* x)
{
    BranchHistory* bh = x->getBranchHistory();

    if (x->getAnc() != NULL) {
        // Only compute mean branch rate if node is NOT the root

        double rate = 0.0;
        int n_events = bh->getNumberOfBranchEvents();

        TraitBranchEvent* ancestralEvent =
            static_cast<TraitBranchEvent*>(bh->getAncestralNodeEvent());

        if (n_events == 0) {

            double t1 = x->getAnc()->getTime();
            double t2 = x->getTime();

            // Times must be relative to event occurrence time:
            t1 -= ancestralEvent->getAbsoluteTime();
            t2 -= ancestralEvent->getAbsoluteTime();

            double zpar = ancestralEvent->getBetaShift();
            double beta0 = ancestralEvent->getBetaInit();

            rate = x->integrateExponentialRateFunction(beta0, zpar, t1, t2);
            rate /= x->getBrlen();

        } else {

            double tcheck = 0.0;
            double t1 = x->getAnc()->getTime();
            double t2 = bh->getEventByIndexPosition(0)->getAbsoluteTime();

            tcheck += (t2 - t1);

            // Times must be relative to initial time of event
            t1 -= ancestralEvent->getAbsoluteTime();
            t2 -= ancestralEvent->getAbsoluteTime();
            double zpar = ancestralEvent->getBetaShift();
            double beta0 = ancestralEvent->getBetaInit();

            rate = x->integrateExponentialRateFunction(beta0, zpar, t1, t2);

            for (int k = 1; k < n_events; k++) {

                t1 = 0.0;
                t2 = bh->getEventByIndexPosition(k)->getAbsoluteTime() -
                     bh->getEventByIndexPosition((k - 1))->getAbsoluteTime();

                TraitBranchEvent* eventAtKMinus1 =
                    static_cast<TraitBranchEvent*>
                        (bh->getEventByIndexPosition(k - 1));

                zpar = eventAtKMinus1->getBetaShift();
                beta0 = eventAtKMinus1->getBetaInit();

                rate += x->integrateExponentialRateFunction(beta0, zpar, t1, t2);

                tcheck += (t2 - t1);

            }

            t1 = 0.0;
            t2 = x->getTime() - bh->getEventByIndexPosition((n_events -
                    1))->getAbsoluteTime();

            TraitBranchEvent* event = static_cast<TraitBranchEvent*>
                (bh->getNodeEvent());

            zpar = event->getBetaShift();
            beta0 = event->getBetaInit();

            rate += x->integrateExponentialRateFunction(beta0, zpar, t1, t2);

            tcheck += (t2 - t1);


            // The overall mean rate across the branch:
            rate /= (x->getBrlen());

            //std::cout << "Rate: " << rate << std::endl;
        }
        x->setMeanBeta(rate);

    } else {
        // Node is root
        x->setMeanBeta((double)0.0);

    }

    TraitBranchEvent* event =
        static_cast<TraitBranchEvent*>(bh->getNodeEvent());

    // compute speciation rate at the focal node:
    double reltime = x->getTime() - event->getAbsoluteTime();

    double init = event->getBetaInit();
    double zz = event->getBetaShift();

    double curBeta = x->getExponentialRate(init, zz, reltime);

#ifdef DEBUG_TIME_VARIABLE

    // Try setting node speciation rates equal to mean rate on descendant branches, to see if the
    //  high-rate trap disappears.

#else

    x->setNodeBeta(curBeta);

#endif

}


/*
 Now setting trait values to be drawn from unifom distribution defined by 2 parental values.

 */
void Tree::initializeTraitValues()
{

    std::cout << "Setting initial trait values at internal nodes" << std::endl;

    // get min & max values:
    double mn = 0;
    double mx = 0;
    bool set = false;

    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getIsTraitFixed()) {
            if (set == false) {
                mn = (*i)->getTraitValue();
                mx = (*i)->getTraitValue();
                set = true;
            } else {
                if ((*i)->getTraitValue() < mn )
                    mn = (*i)->getTraitValue();
                if ((*i)->getTraitValue() > mx)
                    mx = (*i)->getTraitValue();
            }
        }
    }
    recursiveSetTraitValues(root, mn, mx);
}


void Tree::recursiveSetTraitValues(Node* x, double mn, double mx)
{
    if (x->getLfDesc() != NULL && x->getRtDesc() != NULL) {
        recursiveSetTraitValues(x->getLfDesc(), mn, mx);
        recursiveSetTraitValues(x->getRtDesc(), mn, mx);

        // choose random number between two descendants.
        double s1 = x->getLfDesc()->getTraitValue();
        double s2 = x->getRtDesc()->getTraitValue();

        if (s1 < s2) {
            //x->setTraitValue(ranPtr->uniformRv(s1, s2));
            x->setTraitValue((s1 + s2) / (double)2);

        } else if (s1 > s2) {
            //x->setTraitValue(ranPtr->uniformRv(s2, s1));
            x->setTraitValue((s1 + s2) / (double)2);
        } else {
            x->setTraitValue(s1);
        }
    } else if (x->getIsTraitFixed() == false) {
        x->setTraitValue(_random.uniform(mn, mx));
    } else {
        // Trait is fixed. Nothing to do.
    }
}


// Returns pointer to node of mrca of 2 taxa, with names
//  A and B.

Node* Tree::getNodeMRCA(const std::string& A, const std::string& B)
{
    //std::cout << "MRCA of " << A << "\t" << B << std::endl;

    Node* nodeA = NULL;
    Node* nodeB = NULL;
    bool Agood = false;
    bool Bgood = false;

    clearTempNodeArray();
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getName() == A) {
            nodeA = (*i);
            Agood = true;
        }
        if ((*i)->getName() == B) {
            nodeB = (*i);
            Bgood = true;
        }
    }

    if (!Agood | !Bgood) {
        log(Error) << "Invalid nodes " << A << " and " << B
            << " sent to Tree::getNodeMRCA(...)\n";
        //exit(1);
        throw;
    }

    passUpFillTempNodeArray(nodeA);
    //std::cout << _tempNodeSet.size() << std::endl;

    bool isFoundMRCA = false;
    while (!isFoundMRCA) {
        nodeB = nodeB->getAnc();
        //std::cout << nodeB << std::endl;
        if (_tempNodeSet.count(nodeB) > 0) {
            //std::cout << nodeB << " in common" << std::endl;
            break;
        }
    }
    //std::cout << std::endl << std::endl;
    //std::cout << nodeB << std::endl;
    //for (std::set<Node*>::iterator i = _tempNodeSet.begin(); i != _tempNodeSet.end(); i++)
    //  std::cout << "From A: " << (*i) << std::endl;
    clearTempNodeArray();
    return nodeB;
}


void Tree::passUpFillTempNodeArray(Node* x)
{
    if (x != NULL) {
        _tempNodeSet.insert(x);
        passUpFillTempNodeArray(x->getAnc());
    }
}


Node* Tree::getNodeByName(const std::string& A)
{
    Node* x = root;
    int count = 0;
    for (std::vector<Node*>::iterator i = _preOrderNodes.begin();
            i != _preOrderNodes.end(); ++i) {
        if ((*i)->getName() == A) {
            count++;
            x = (*i);
        }
    }
    if (count == 0) {
        std::cout << "Invalid node name: name not found in Tree:: getNodeByName" << std::endl;
        exit(0);
    } else if (count > 1) {
        std::cout << "Duplicate node names found in Tree:: getNodeByName" << std::endl;
        exit(0);
    } //else {

    //}

    return x;
}


double Tree::maxRootToTipLength()
{
    std::vector<double> lengths = terminalPathLengthsToRoot();
    return *std::max_element(lengths.begin(), lengths.end());
}


bool Tree::isUltrametric()
{
    std::vector<double> terminalPathLengths = terminalPathLengthsToRoot();
    return Stat::variance(terminalPathLengths) < ULTRAMETRIC_TOLERANCE;
}


std::vector<double> Tree::terminalPathLengthsToRoot()
{
    std::vector<double> pathLengths;
    storeTerminalPathLengthsToRootRecurse(getRoot(), pathLengths);
    return pathLengths;
}


void Tree::storeTerminalPathLengthsToRootRecurse
    (Node* node, std::vector<double>& pathLengths)
{
    Node* leftNode = node->getLfDesc();
    Node* rightNode = node->getRtDesc();

    if (leftNode == NULL && rightNode == NULL) {
        pathLengths.push_back(node->pathLengthToRoot());
        return;
    }

    if (leftNode != NULL) {
        storeTerminalPathLengthsToRootRecurse(leftNode, pathLengths);
    }

    if (rightNode != NULL) {
        storeTerminalPathLengthsToRootRecurse(rightNode, pathLengths);
    }
}


void Tree::assertTreeRootBranchLengthIsZero()
{
    Node* root = getRoot();

    if (root->getBrlen() != 0.0) {
        log(Warning) << "Root has non-zero branch length. "
            << "Automatically setting it to zero.\n"
            << "Please make sure this is what you really want.\n\n";
        root->setBrlen(0.0);
    }
}


void Tree::assertTreeIsBifurcating()
{
    assertTreeIsBifurcatingRecurse(getRoot());
}


void Tree::assertTreeIsBifurcatingRecurse(Node* node)
{
    Node* left = node->getLfDesc();
    Node* right = node->getRtDesc();

    if ((left == NULL && right != NULL) ||
        (left != NULL && right == NULL)) {
        log(Error) << "Tree is not bifurcating.\n"
            << "The node with branch length " << node->getBrlen() << " "
            << "has only one child node.\n";
        std::exit(1);
    }

    if (left != NULL && right != NULL) {
        assertTreeIsBifurcatingRecurse(left);
        assertTreeIsBifurcatingRecurse(right);
    }
}


void Tree::assertBranchLengthsArePositive()
{
    Node* root = getRoot();

    if (root->getBrlen() < 0.0) {
        log(Error) << "Branch length of root node is negative.\n";
        std::exit(1);
    }

    assertBranchLengthsArePositiveRecurse(root->getLfDesc());
    assertBranchLengthsArePositiveRecurse(root->getRtDesc());
}


void Tree::assertBranchLengthsArePositiveRecurse(Node* node)
{
    if (node == NULL) {
        return;
    }

    if (node->getBrlen() <= 0.0) {
        log(Error) << "At least one branch length is non-positive.\n";
        std::exit(1);
    }

    assertBranchLengthsArePositiveRecurse(node->getLfDesc());
    assertBranchLengthsArePositiveRecurse(node->getRtDesc());
}


void Tree::assertTreeIsUltrametric()
{
    if (!isUltrametric()) {
        log(Error) << "Tree is not ultrametric.\n";
        std::exit(1);
    }
}


void Tree::assertTipsHaveUniqueNames()
{
    std::vector<std::string> names = terminalNames();

    // Sort names
    std::sort(names.begin(), names.end());

    // Find consecutive duplicates
    std::vector<std::string>::const_iterator it;
    it = std::adjacent_find(names.begin(), names.end());

    if (it != names.end()) {
        log(Error) << "Tree contains tips with the same name.\n";
        std::exit(1);
    }
}


std::vector<std::string> Tree::terminalNames()
{
    std::vector<std::string> names;
    storeTerminalNamesRecurse(getRoot(), names);
    return names;
}


void Tree::storeTerminalNamesRecurse
    (Node* node, std::vector<std::string>& names)
{
    Node* leftNode = node->getLfDesc();
    Node* rightNode = node->getRtDesc();

    if (leftNode == NULL && rightNode == NULL) {
        names.push_back(node->getName());
        return;
    }

    if (leftNode != NULL) {
        storeTerminalNamesRecurse(leftNode, names);
    }

    if (rightNode != NULL) {
        storeTerminalNamesRecurse(rightNode, names);
    }
}


// TODO: Use this function in the functions above
std::vector<Node*> Tree::terminalNodes()
{
    std::vector<Node*> nodes;
    storeTerminalNodesRecurse(getRoot(), nodes);
    return nodes;
}


void Tree::storeTerminalNodesRecurse(Node* node, std::vector<Node*>& nodes)
{
    Node* leftNode = node->getLfDesc();
    Node* rightNode = node->getRtDesc();

    if (leftNode == NULL && rightNode == NULL) {
        nodes.push_back(node);
        return;
    }

    if (leftNode != NULL) {
        storeTerminalNodesRecurse(leftNode, nodes);
    }

    if (rightNode != NULL) {
        storeTerminalNodesRecurse(rightNode, nodes);
    }
}


std::vector<double> Tree::traitValues()
{
    const std::vector<Node*>& nodes = terminalNodes();
    std::vector<double> values;

    std::vector<Node*>::const_iterator it;
    for (it = nodes.begin(); it != nodes.end(); ++it) {
        values.push_back((*it)->getTraitValue());
    }

    return values;
}
