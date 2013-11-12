#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "Tree.h"
#include "Node.h"

#include "BranchHistory.h"
#include "TraitBranchHistory.h"
#include "MbRandom.h"
//#include "phenotype.h"


//#define DEBUG_TIME_VARIABLE

Tree::Tree(void)
{

}

Tree::~Tree(void)
{
    //cout << "calling destructor" << endl;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        delete (*i);

}

Tree::Tree(string fname, MbRandom* rnptr)
{
    //cout << "in constructor..." << endl;

    ranPtr = rnptr;

    string treestring;

    ifstream treefile(fname.c_str());

    //treefile.open(fname.c_str());
    cout << "Reading tree from file <" << fname << ">" << endl;

    if (!treefile.good()) {
        cout << "Invalid filename for phylogenetic tree\n" << endl;
        throw;

    }


    treefile >> treestring;
    treefile.close();

    //cout << "tree size: " << treestring.size() << endl;

    setTaxonCountFromNewickString(treestring);

    buildTreeFromNewickString(treestring);

    getDownPassSeq();

    // All of this below is output and can be deleted:

    // Output stuff here
    cout << "1 tree read with " << getNumberTips() << " taxa" << endl;


    // counting tips for trial...
    int sum = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        sum += (int)(*i)->getIsTip();
    //cout << sum << endl;

    // initialize treelength:
    treelength = 0.0;
    _totalMapLength = 0.0;

    // Set node times (with 0 at root):
    setNodeTimes(root);
    setAge();
    setBranchingTimes(root);


    //Node* tmp = root->getRtDesc();
    //cout << "tmptime: " << tmp->getMapStart() << "\t" << tmp->getMapEnd() << endl;

    //cout << "roottime: " << root->getMapStart() << "\t" << root->getMapEnd() << endl;

    //cout << root->getBrlen() << " root brlen" << endl;

    // Need to set treelength:
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        treelength += (*i)->getBrlen();

    //cout << "treelength: " << treelength << endl << endl;

    // Setting internal node set:
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ( (*i)->getLfDesc() != NULL && (*i)->getRtDesc() != NULL )
            internalNodeSet.insert((*i));
    }

    // Set tip counts for each node.
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        int dcount = getDescTipCount((*i));
        (*i)->setTipDescCount(dcount);
        // cout << (*i)->getLfDesc() << "\t" << (*i)->getRtDesc() << "\t";
        // cout << (*i)->getTipDescCount() << "\tCanHold: " << (*i)->getCanHoldEvent() << "\t" << (*i) << endl;
    }


    //cout << "Number of internal nodes: " << internalNodeSet.size() << endl << endl;
    int ct = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getCanHoldEvent())
            ct++;
    }
    cout << "Tree ctor: event nodes: " << ct << endl;

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
        if (p->getRtDesc()->getCanHoldEvent())
            setTreeMap(p->getRtDesc());
    }
    if (p->getLfDesc() != NULL) {
        if (p->getLfDesc()->getCanHoldEvent())
            setTreeMap(p->getLfDesc());
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
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {

        if (x >= (*i)->getMapStart()  && x < (*i)->getMapEnd()) {


            // Very serious bug, fixed 9.15.2012
            //double delta = (*i)->getMapEnd() - x;
            //abstime = (*i)->getTime() - delta;
            //done = true;
            // code above should be inverting event positions on branch.


            double delta = x - (*i)->getMapStart(); // difference in times...
            abstime = (*i)->getTime() - delta;
            done = true;
        }

    }
    if (done == false) {
        cout << "could not find abs time from map time \n";
        cout << "Tree::getAbsoluteTimeFromMapTime() " << endl;
        throw;
    }


    return abstime;
}

void Tree::printCanHoldEventByNode(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        cout << (*i) << "\tCan hold: " << (*i)->getCanHoldEvent() << "\tTips: " <<
             (*i)->getTipDescCount() << endl;

}

// Get number of descendant nodes from a given node
int Tree::getDescNodeCount(Node* p)
{
    double count = 0;
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
        if ( x > (*i)->getMapStart() && x <= (*i)->getMapEnd())
            y = (*i);
    }
    if (y == NULL) {
        cout << "error: unmapped event\n" << endl;
        cout << "position: " << x << endl;
    }
    //cout << "y in map: " << y << endl;
    return y;
}

void Tree::printNodeMap(void)
{
    for (std::set<Node*>::iterator i = mappableNodes.begin();
            i != mappableNodes.end(); i++)
        cout << (*i) << "\t" << (*i)->getAnc() << "\t" << (*i)->getMapStart() << "\t"
             << (*i)->getMapEnd() << endl;



}


void Tree::getDownPassSeq(void)
{

    passDown(root);

}


void Tree::passDown(Node* p)
{
    if (p != NULL) {
        passDown(p->getLfDesc());
        passDown(p->getRtDesc());
        downPassSeq.push_back(p);
    }


}


// This requires that p be an internal node to begin with!
void Tree::setTempInternalNodeArray(Node* p)
{

    if (p->getRtDesc() == NULL && p->getLfDesc() == NULL) {
        cout << "Problem: sent terminal node to setTempInternalNodeArray" << endl;
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

void Tree::clearTempNodeArray(void)
{
    //cout << "Size before: " << _tempNodeSet.size() << endl;

    for (std::set<Node*>::iterator i = _tempNodeSet.begin();
            i != _tempNodeSet.end(); i++)
        _tempNodeSet.erase(i);
    //cout << "Size after: " << _tempNodeSet.size() << endl;
}

Node* Tree::getRandomNodeFromTempArray(void)
{

    int chosen = ranPtr->sampleInteger(0, (_tempNodeSet.size() - 1));
    int myit = 0;
    Node* xnode = (*_tempNodeSet.begin());

    for (std::set<Node*>::iterator i = _tempNodeSet.begin();
            i != _tempNodeSet.end(); i++) {
        if (myit == chosen)
            xnode = (*i);
        myit++;

    }

    return xnode;

}


void Tree::setBranchLengths(void)
{

    // Set brlens:
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i) == root) {
            double timeX = (*i)->getTime() - _startTime;
            (*i)->setBrlen( timeX );
        } else {
            double timeX = (*i)->getTime() - ( (*i)->getAnc() )->getTime();
            (*i)->setBrlen( timeX );
        }

    }
}


string Tree::getNewick(void)
{

    stringstream ss;

    writeTree(root, ss);
    string newick = ss.str();
    newick.append(";");
    return newick;

}


void Tree::writeTree(Node* p, stringstream& ss)
{

    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "")
            ss << p->getIndex() << ":" << p->getBrlen();
        else
            ss << p->getName() << ":" << p->getBrlen();
    } else {
        ss << "(";
        writeTree(p->getLfDesc(), ss);
        ss << ",";
        writeTree(p->getRtDesc(), ss);
        ss << "):" << p->getBrlen();
    }

}


int Tree::getNumberTips(void)
{

    int count = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ( (*i)->getDescCount() == 0 )
            count++;
    }
    return count;
}

void Tree::setTipStatus(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if (((*i)->getLfDesc() == NULL) && ((*i)->getRtDesc() == NULL ))
            (*i)->setIsTip(true);
        else
            (*i)->setIsTip(false);
    }

}

int Tree::getNumberExtantTips(void)
{
    setExtantStatus();
    int ntaxa = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if (((*i)->getLfDesc() == NULL) && ((*i)->getExtantStatus() == 1) )
            ntaxa++;
    }
    return ntaxa;
}



/* setExtantStatus
        value = 1 if node is a leaf AND if
        node is extant
        OR if node is internal

        extinct leaves get value = 0

 */


void Tree::setExtantStatus(void)
{
    double tx = root->getTime();
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getTime() > tx)
            tx = (*i)->getTime();
    }

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getLfDesc() != NULL && (*i)->getRtDesc() != NULL)
            (*i)->setExtantStatus(1);
        else {
            if ((*i)->getTime() > (tx - 0.000001) )
                (*i)->setExtantStatus(1);
            else
                (*i)->setExtantStatus(0);
        }

        //cout << "Node\t" << (*i)->getIndex() << "\t" << (*i)->getExtantStatus() << endl;
    }

}


/*
 This function iterates over all nodes.
 If the node is a tip, and ALIVE, it gets flagged with 1. Otherwise 0.
 This can be distinguished from setExtantStatus,
    which flags status of internal nodes as well.
 */
void Tree::setIsLivingTipStatus(void)
{

    double TOL = 0.000001;

    // Get maximum age.
    double tx = root->getTime();
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getTime() > tx)
            tx = (*i)->getTime();
    }

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL) {
            double test = abs((*i)->getTime() - tx );
            if (test < TOL)
                (*i)->setIsLivingTip(true);
            else
                (*i)->setIsLivingTip(false);
        } else
            (*i)->setIsLivingTip(false);
    }

}

void Tree::writeNodeData(void )
{

    int wsize = 15;

    ofstream myfile;
    myfile.open("NodeData.txt");
    myfile << "index" << setw(wsize) << "LDindex" << setw(wsize);
    myfile << "RDindex" << setw(wsize) << "Node" << setw(wsize);
    myfile << "LDesc" << setw(wsize) << "RDesc" << setw(wsize);
    myfile << "anc" << setw(wsize) << "time" << endl;


    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {

        myfile << (*i)->getIndex() << setw(wsize);

        if ( (*i)->getLfDesc() != 0 ) {
            myfile << (*i)->getLfDesc()->getIndex() << "\t";
            myfile << (*i)->getRtDesc()->getIndex() << "\t";
        } else {
            myfile << -1 << "\t";
            myfile << -1 << "\t";
        }

        myfile << (*i) << setw(wsize);
        myfile << (*i)->getRtDesc() << setw(wsize);
        myfile << (*i)->getLfDesc() << setw(wsize);
        myfile << (*i)->getAnc() << setw(wsize);
        myfile << (*i)->getTime() << setw(wsize);
        myfile << endl;

    }

    myfile.close();

}

void Tree::rebuildTreeNodeSet(void)
{
    nodes.clear();
    recursivelyAddNodesToSet(root);

}

void Tree::recursivelyAddNodesToSet(Node* p)
{
    if (p != NULL) {
        nodes.insert(p);
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



void Tree::pruneExtinctTaxa(void)
{
    //assume most recent time is extant

    setExtantStatus();

    fixExtinct(root);


    //cout << nodes.size() << endl;
    rebuildTreeNodeSet();
    //cout << nodes.size() << endl;

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

    if (p->getDescCount() == 0)
        cout << "terminal node: should never get here\n" << endl;
    // separate recursion from modification:
    // Recursion step:
    if (p->getLfDesc()->getDescCount() == 2)
        fixExtinct(p->getLfDesc());
    if (p->getRtDesc()->getDescCount() == 2)
        fixExtinct(p->getRtDesc());

    // modification step:

    if (p->getLfDesc()->getExtantStatus() == 0 &&
            p->getRtDesc()->getExtantStatus() == 1) {
        //cout << "here" << endl;
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

        } else
            cout << "problem: p state invalid" << endl;


    } else if (p->getLfDesc()->getExtantStatus() == 1 &&
               p->getRtDesc()->getExtantStatus() == 0) {

        if (p == root) {
            //cout << "p is root" << endl;
            //cout << root << "\t" << p << "\t" << endl;
            root = p->getLfDesc();
            //cout << p->getLfDesc() << "\t" << root << endl;
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
        } else
            cout << "invalid p state" << endl;

    } else if (p->getLfDesc()->getExtantStatus() == 0 &&
               p->getRtDesc()->getExtantStatus() == 0) {
        if (p == root)
            cout << "Problem: tree extinct" << endl;
        else {
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
    } else
        cout << "problem: invalid condition" << endl;
}





/* deleteExtinctNodes
 This function does:
    1. find extinct node
    2. set pointer pointing to node to NULL
    3. delete node

 */

void Tree::deleteExtinctNodes(void)
{

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        // if node is leaf:
        cout << (*i)->getIndex() << endl;
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL) {
            // if nodes is extinct
            if ((*i)->getExtantStatus() == 0) {
                if ((*i)->getAnc()->getLfDesc() == (*i)) {
                    (*i)->getAnc()->setLfDesc(NULL);
                    //delete (*i);
                } else if ((*i)->getAnc()->getRtDesc() == (*i) ) {
                    (*i)->getAnc()->setRtDesc(NULL);
                    //delete (*i);
                } else
                    cout << "problem in deleteExtinctNodes()" << endl;

            }
        }


    }

}


void Tree::setTaxonCountFromNewickString(string ts)
{
    int count = 0;
    for (string::size_type i = 0; i < ts.size(); i++) {
        char c = ts[i];
        if (c == ',')
            count++;
    }
    _ntaxa = count + 1;
}


/*
 Rewrite this one to build node array as we go along
 through string
 */


/*I think this works...*/

void Tree::buildTreeFromNewickString(string ts)
{

    //cout << "in build tree..." << endl;

    bool readingBL = false;
    Node* p = NULL;

    //int nextInterNode = _ntaxa;
    //int taxCounter = 0;
    //int nodecounter = 0;

    //std::set<Node*>::iterator NodeIterator = nodes.begin();

    for (string::size_type i = 0; i < ts.size(); i++) {

        char c = ts[i];
        //cout << c << endl;
        if (c == '(') {
            //q = &nodes[nextInterNode++];

            Node* q = new Node;
            nodes.insert(q);

            //q = *NodeIterator++;
            //cout << ++nodecounter << '\t' << p << '\t' << q << endl;

            if (p == NULL) {
                p = q;
                root = p;
            } else {
                //cout << p->getLfDesc() << "\t" << p->getRtDesc() << endl;
                q->setAnc(p);
                if (p->getLfDesc() == NULL)
                    p->setLfDesc(q);
                else if (p->getRtDesc() == NULL)
                    p->setRtDesc(q);
                else {
                    cerr << "ERROR: tree string";
                    exit(1);
                }
            }
            p = q;
            readingBL = false;
        } else if (c == ')') {
            if (p->getAnc() == NULL) {
                cerr << "ERROR: tree string";
                exit(1);
            } else p = p->getAnc();
            readingBL = false;
        } else if (c == ',') {
            if (p->getAnc() == NULL) {
                cerr << "ERROR: tree string";
                exit(1);
            } else p = p->getAnc();
            readingBL = false;
        } else if (c == ':')
            readingBL = true;
        else if (c == ';') {
            // done with tree
        } else {
            string s = "";
            while (isValidChar(ts[i]))
                s += ts[i++];
            i--;
            if (readingBL == false) {
                // set tip name

                //q = &nodes[taxCounter];
                //q = *NodeIterator++;
                //cout << ++nodecounter << endl;
                Node* q = new Node();
                nodes.insert(q);

                if (p == NULL) {
                    cerr << "ERROR: Problem adding a tip to the tree" << endl;
                    exit(1);
                } else {
                    q->setAnc(p);
                    if (p->getLfDesc() == NULL)
                        p->setLfDesc(q);
                    else if (p->getRtDesc() == NULL)
                        p->setRtDesc(q);
                    else {
                        cerr << "ERROR: Problem adding a tip to the tree" << endl;
                        exit(1);
                    }
                }
                p = q;
                p->setName(s);
                p->setIsTip(true);
                //taxCounter++;
            } else {
                // read in bl
                double v = 0.0;
                istringstream buf(s);
                buf >> v;
                p->setBrlen(v);
                //p->setSimmedBrLen(v);
            }
        }
    }

    setStartTime(0);
    setNodeTimes(root);
    setAge();


}

/* sets time of each node*/

void Tree::setNodeTimes(Node* p)
{
    if (p == root)
        p->setTime(0);
    else {
        double x = p->getBrlen() + p->getAnc()->getTime();
        p->setTime(x);
    }

    if (p->getLfDesc() != NULL && p->getRtDesc() != NULL) {
        setNodeTimes(p->getLfDesc());
        setNodeTimes(p->getRtDesc());
    }

}


bool Tree::isValidChar(char x)
{
    if (x == ')' || (x == '(') || (x == ':') || (x == ',') || (x == ';'))
        return false;
    else
        return true;
}


void Tree::setAge(void)
{
    // age here defined as MAX node time (node times start with 0 at root)
    double mx = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getTime() == 0)
            setNodeTimes(root);
        if ((*i)->getTime() > mx)
            mx = (*i)->getTime();
    }
    _age = mx;
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

vector<double>  Tree::getBranchingTimes(void)
{

    double TOL = 0.0000001;

    vector<double> btimes;

    //cout << "age :" << _age << endl;

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        double tmp = _age - (*i)->getTime();
        //cout << "Time: " << (*i)->getTime() << "\t" << tmp << endl;

        if (tmp > TOL)
            btimes.push_back(tmp);
    }
    sort(btimes.begin(), btimes.end() );
    reverse(btimes.begin(), btimes.end());
    return btimes;
}


void Tree::writeMeanBranchTraitRateTree(Node* p, stringstream& ss)
{

    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "")
            ss << p->getIndex() << ":" << p->getMeanBeta();
        else
            ss << p->getName() << ":" << p->getMeanBeta();
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


// Update both mean speciation rates on branch in addition to node speciation rate.
void Tree::setMeanBranchSpeciation(void)
{

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        (*i)->computeNodeBranchSpeciationParams();
}

// this will update extinction for both branch mean rates as well as
// node-associated rates.

void Tree::setMeanBranchExtinction(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        (*i)->computeNodeBranchExtinctionParams();
}

void Tree::echoMeanBranchRates(void)
{

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        Node* x = (*i);
        cout << x << "\t" << x->getMeanSpeciationRate() << "\t" <<
             x->getMeanExtinctionRate() << "\t" << x->getNodeLambda() << endl;

    }


}

void Tree::echoMeanBranchTraitRates(void)
{

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        Node* x = (*i);
        cout << x->getName() << "\t" << x->getMeanBeta() << endl;

    }


}



void Tree::writeBranchSpeciationRatesToFile(string fname, bool append)
{


    stringstream outdata;

    writeMeanBranchSpeciationTree(root, outdata);

    outdata << ";";

    ofstream outStream;
    if (append == true)
        outStream.open(fname.c_str(), ofstream::app);
    else
        outStream.open(fname.c_str(), ofstream::trunc);
    outStream << outdata.str() << endl;
    outStream.close();

}


void Tree::writeBranchExtinctionRatesToFile(string fname, bool append)
{


    stringstream outdata;

    writeMeanBranchExtinctionTree(root, outdata);

    outdata << ";";

    ofstream outStream;
    if (append == true)
        outStream.open(fname.c_str(), ofstream::app);
    else
        outStream.open(fname.c_str(), ofstream::trunc);
    outStream << outdata.str() << endl;
    outStream.close();

}


void Tree::writeMeanBranchSpeciationTree(Node* p, stringstream& ss)
{

    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "")
            ss << p->getIndex() << ":" << p->getMeanSpeciationRate();
        else
            ss << p->getName() << ":" << p->getMeanSpeciationRate();
    } else {
        ss << "(";
        writeMeanBranchSpeciationTree(p->getLfDesc(), ss);
        ss << ",";
        writeMeanBranchSpeciationTree(p->getRtDesc(), ss);
        ss << "):" << p->getMeanSpeciationRate();
    }


}

void Tree::writeNodeSpeciationTree(Node* p, stringstream& ss)
{

    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "")
            ss << p->getIndex() << ":" << p->getNodeLambda();
        else
            ss << p->getName() << ":" << p->getNodeLambda();
    } else {
        ss << "(";
        writeNodeSpeciationTree(p->getLfDesc(), ss);
        ss << ",";
        writeNodeSpeciationTree(p->getRtDesc(), ss);
        ss << "):" << p->getNodeLambda();
    }

}



void Tree::writeMeanBranchExtinctionTree(Node* p, stringstream& ss)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "")
            ss << p->getIndex() << ":" << p->getMeanExtinctionRate();
        else
            ss << p->getName() << ":" << p->getMeanExtinctionRate();
    } else {
        ss << "(";
        writeMeanBranchExtinctionTree(p->getLfDesc(), ss);
        ss << ",";
        writeMeanBranchExtinctionTree(p->getRtDesc(), ss);
        ss << "):" << p->getMeanExtinctionRate();
    }

}

void Tree::writeMeanBranchNetDivRateTree(Node* p, stringstream& ss)
{

    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "")
            ss << p->getIndex() << ":" << (p->getMeanSpeciationRate() -
                                           p->getMeanExtinctionRate());
        else
            ss << p->getName() << ":" << (p->getMeanSpeciationRate() -
                                          p->getMeanExtinctionRate());
    } else {
        ss << "(";
        writeMeanBranchNetDivRateTree(p->getLfDesc(), ss);
        ss << ",";
        writeMeanBranchNetDivRateTree(p->getRtDesc(), ss);
        ss << "):" << (p->getMeanSpeciationRate() - p->getMeanExtinctionRate());
    }




}

void Tree::writeBranchPhenotypes(Node* p, stringstream& ss)
{
    if (p->getLfDesc() == NULL && p-> getRtDesc() == NULL) {
        if (p->getName() == "")
            ss << p->getIndex() << ":" << p->getTraitValue();
        else
            ss << p->getName() << ":" << p->getTraitValue();
    } else {
        ss << "(";
        writeBranchPhenotypes(p->getLfDesc(), ss);
        ss << ",";
        writeBranchPhenotypes(p->getRtDesc(), ss);
        ss << "):" << p->getTraitValue();
    }


}


/*
Tree::getPhenotypes
Read file. First column = species name exactly as matching in phylogeny.
 Second column, tab-delimited: phenotype value on appropriate scale
    (e.g., already log-transformed).

 */

void Tree::getPhenotypes(string fname)
{

    ifstream infile(fname.c_str());
    cout << "Reading phenotypes from file <<" << fname.c_str() << ">>" << endl;
    vector<string> stringvec;

    vector<string> spnames;
    vector<double> traits;

    //treefile.open(fname.c_str());

    if (!infile.good())
        cout << "Bad Filename" << endl;

    while (infile) {
        string tempstring;
        getline(infile, tempstring, '\t');

        //cout << tempstring << "\n" << endl;

        spnames.push_back(tempstring);

        getline(infile, tempstring, '\n');

        traits.push_back(atof(tempstring.c_str()));

        // this OK?
        if (infile.peek() == EOF)
            break;

    }


    infile.close();

    cout << "Read a total of " << traits.size() << " species w trait data" << endl;

    // iterate over nodes...
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL ) {

            for (vector<string>::size_type k = 0; k < spnames.size(); k++) {
                if ((*i)->getName() == spnames[k]) {
                    (*i)->setTraitValue(traits[k]);
                    (*i)->setIsTraitFixed(true);

                }

            }
            if ((*i)->getIsTraitFixed() == false)
                cout << "error - failed to set a terminal state\n" << endl;

        } else {
            (*i)->setTraitValue(0);
            (*i)->setIsTraitFixed(false);
        }



    }


}



void Tree::getPhenotypesMissingLatent(string fname)
{

    ifstream infile(fname.c_str());
    cout << "Reading phenotypes from file <<" << fname.c_str() << ">>" << endl;
    vector<string> stringvec;

    vector<string> spnames;
    vector<double> traits;

    //treefile.open(fname.c_str());

    if (!infile.good())
        cout << "Bad Filename" << endl;

    while (infile) {
        string tempstring;
        getline(infile, tempstring, '\t');

        //cout << tempstring << "\n" << endl;

        spnames.push_back(tempstring);

        getline(infile, tempstring, '\n');

        traits.push_back(atof(tempstring.c_str()));

        // this OK?
        if (infile.peek() == EOF)
            break;

    }

    int missingTerminalCount = 0;

    infile.close();

    cout << "Read a total of " << traits.size() << " species w trait data" << endl;

    // iterate over nodes...
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL ) {

            for (vector<string>::size_type k = 0; k < spnames.size(); k++) {
                if ((*i)->getName() == spnames[k]) {
                    (*i)->setTraitValue(traits[k]);
                    (*i)->setIsTraitFixed(true);
                }

            }
            if ((*i)->getIsTraitFixed() == false)
                missingTerminalCount++;

        } else {
            (*i)->setTraitValue(0);
            (*i)->setIsTraitFixed(false);
        }



    }

    cout << "Missing data for < " << missingTerminalCount << " > species." << endl;
    cout << "These will be treated as latent variables in analysis" << endl << endl;

    int count2 = 0;
    int count3 = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getIsTraitFixed())
            count2++;
        else
            count3++;
    }
    cout << "count of FIXED nodes in getPhenotypesMissingLatent: " << count2 <<
         endl;
    cout << "count of VARIABLE nodes in getPhenotypesMissingLatent: " << count3 <<
         endl << endl;
}




void Tree::printTraitValues(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        cout << (*i) << "\t" << (*i)->getIsTraitFixed() << "\t" <<
             (*i)->getTraitValue() << endl;

}

// recursive function to generate trait values...
void Tree::generateTraitsAllNodesBM(Node* xnode, double varx)
{


    if (xnode->getLfDesc() != NULL && xnode->getRtDesc() != NULL) {
        // do left:
        double vx = xnode->getLfDesc()->getBrlen() * varx;
        double newTrait = xnode->getTraitValue() + ranPtr->normalRv((double)0.0,
                          sqrt(vx));
        xnode->getLfDesc()->setTraitValue(newTrait);
        generateTraitsAllNodesBM(xnode->getLfDesc(), varx);


        vx = xnode->getRtDesc()->getBrlen() * varx;
        newTrait = xnode->getTraitValue() + ranPtr->normalRv((double)0.0, sqrt(vx));
        xnode->getRtDesc()->setTraitValue(newTrait);
        generateTraitsAllNodesBM(xnode->getRtDesc(), varx);


    } else if (xnode->getLfDesc() == NULL && xnode->getRtDesc() == NULL)
        xnode->setIsTraitFixed(true);
    else {
        cout << "error in Tree::generateTraitsAllNodesBM" << endl;
        throw;
    }

    //cout << xnode->getTraitValue() << endl;
}






/*
    chooseInternalNodeAtRandom()
        have checked this to make sure distribution of sampled nodes is uniform
        appears to be fine.

 */

Node* Tree::chooseInternalNodeAtRandom(void)
{


    int snode = ranPtr->sampleInteger((int)1, (int)internalNodeSet.size());

    std::set<Node*>::iterator myIt = internalNodeSet.begin();

    for (int i = 1; i < snode; i++ )
        myIt++;
    return  (*myIt);

}



/*
    This sets all internal nodes equal to the estimated sampling fraction for that particular node.
    Ei values for all tip nodes get set;
    Di values for all tip nodes get set.
    Etip gets set for all nodes for global sampling probability

*/
void Tree::initializeSpeciationExtinctionModel(double sampFrac)
{

    double  speciationInit = sampFrac;
    double  extinctionInit = (double)1 - sampFrac;

    for (std::set<Node*>::iterator myIt = nodes.begin(); myIt != nodes.end();
            myIt++) {
        if ((*myIt)->getLfDesc() == NULL && (*myIt)->getRtDesc() == NULL) {
            (*myIt)->setDinit(speciationInit);
            (*myIt)->setEinit(extinctionInit);
        }
        (*myIt)->setEtip(extinctionInit); // Set
    }

    cout << "Speciation/Extinction initial conditions set (global)" << endl;

}

void Tree::initializeSpeciationExtinctionModel(string fname)
{
    // Part 1. Read from file

    ifstream infile(fname.c_str());
    cout << "Reading sampling fractions from file <<" << fname.c_str() << ">>" <<
         endl;
    vector<string> stringvec;

    vector<string> spnames;
    vector<string> spfamilies;
    vector<double> sfracs;

    if (!infile.good())
        cout << "Bad Filename" << endl;

    string tempstring;

    // First number in file is sampling probability for "backbone" of the tree.

    getline(infile, tempstring, '\n');
    double backboneSampProb = atof(tempstring.c_str());
    double backboneInitial = 1 - backboneSampProb;

    while (infile) {
        //string tempstring;
        getline(infile, tempstring, '\t');

        spnames.push_back(tempstring);

        getline(infile, tempstring, '\t');

        spfamilies.push_back(tempstring);

        getline(infile, tempstring, '\n');

        sfracs.push_back(atof(tempstring.c_str()));

        if (infile.peek() == EOF)
            break;

    }

    infile.close();

    cout << "Read a total of " << sfracs.size() << " initial vals" << endl;

    int counter = 0;
    for (std::vector<Node*>::iterator i = downPassSeq.begin();
            i != downPassSeq.end(); i++) {

        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL ) {
            for (vector<string>::size_type k = 0; k < spnames.size(); k++) {
                if ((*i)->getName() == spnames[k]) {

                    double Einit = (double)1 - sfracs[k];
                    double Dinit = sfracs[k];

                    (*i)->setEinit(Einit);
                    (*i)->setEtip(Einit);
                    (*i)->setDinit(Dinit);
                    (*i)->setCladeName(spfamilies[k]);
                    counter++;
                }
                //cout << spfamilies[k] << endl;
            }
            if ((*i)->getEinit() == -1)
                cout << ((*i)->getName()) << endl;
        } else {
            // Node is internal
            if ((*i)->getLfDesc()->getCladeName() == (*i)->getRtDesc()->getCladeName()) {
                // node *i belongs to same clade and inherits their sampling probability:
                double sprob = (*i)->getLfDesc()->getEtip();
                (*i)->setEtip(sprob);
                (*i)->setCladeName((*i)->getLfDesc()->getCladeName());
            } else {
                string cname = "backbone";
                (*i)->setCladeName(cname);
                (*i)->setEtip(backboneInitial);
            }


        }
        //cout << "here?" << endl;
        if ((*i)->getCladeName() == "") {
            cout << "unset clade names \n" << endl;
            throw;
        }

    }

    /*

    // Part 2. Do the initialization for terminals.
    int counter = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++){
        if ((*i)->getLfDesc() == NULL && (*i)->getRtDesc() == NULL ){
            for (int k = 0; k < spnames.size(); k++){
                if ((*i)->getName() == spnames[k]){

                    double Einit = (double)1 - sfracs[k];
                    double Dinit = sfracs[k];

                    (*i)->setEinit(Einit);
                    (*i)->setDinit(Dinit);
                    counter++;
                }

            }
            if((*i)->getEinit() == -1){
                cout << ((*i)->getName()) << endl;
            }
        }
    }
    */

    /*
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++){
        Node * x = (*i);
        cout << x->getLfDesc() << "\t" << x->getEtip() << "\t" << x->getCladeName() << endl;
        //cout << x->getName() << "\t" << x->getEinit() << "\t" << x->getEtip() << "\t" << x->getCladeName() << "\t";
        cout << endl;
    }
    */

    int tcount = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        Node* x = (*i);
        if (x->getEtip() < 0)
            tcount++;

        //cout << x->getLfDesc() << "\t" << x->getEtip() << "\t" << x->getCladeName() << endl;
        //cout << x->getName() << "\t" << x->getEinit() << "\t" << x->getEtip() << "\t" << x->getCladeName() << "\t";
        //cout << endl;
    }


    cout << "Set a total of < " << counter << " > tips nodes for Ei & Di" << endl;
    cout << "Failed to set < " << tcount << " > internal node eTip values" << endl;

}



void Tree::printInitialSpeciationExtinctionRates(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        //if ((*i)->getLfDesc() == NULL && ((*i)->getRtDesc() == NULL))
        cout << (*i) << "\t" << (*i)->getDinit() << "\t" << (*i)->getEinit()  << "\t"
             << (*i)->getBrlen() << endl;
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

    for (std::set<Node*>::iterator myIt = nodes.begin(); myIt != nodes.end();
            myIt++) {
        int dcount = countDescendantsWithValidTraitData((*myIt));
        if (dcount >= ndesc) {
            (*myIt)->setCanHoldEvent(true);
            mappableNodes.insert((*myIt));
        }
    }

    cout << "Number of mappable nodes: < " << mappableNodes.size() << " >" << endl;

}

int Tree::getDescTipCount(Node* p)
{
    int count = 0;
    if (p->getLfDesc() == NULL & p->getRtDesc() == NULL)
        count++;
    else {
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
    if (p->getLfDesc() == NULL & p->getRtDesc() == NULL) {
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

void Tree::setAllNodesCanHoldEvent(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        (*i)->setCanHoldEvent(true);
        mappableNodes.insert((*i));
    }

}

void Tree::setCanNodeHoldEventByDescCount(int x)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
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
    cout << "Map length: " << getTotalMapLength() << endl;
    cout << "Total mappable nodes: " << mappableNodes.size() << endl;

}


void Tree::printTraitRange(void)
{

    double mx = root->getTraitValue();
    double mn = root->getTraitValue();

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        double ctrait = (*i)->getTraitValue();
        if (ctrait < mn)
            mn = ctrait;
        if (ctrait > mx)
            mx = ctrait;


    }

    cout << "Min trait value < " << mn << " >\tMax value < " << mx << " >" << endl;

}

void Tree::loadPreviousNodeStates(Tree* ostree)
{

    cout << "ostree read into loadpreviousstates" << endl;


    for (vector<Node*>::size_type i = 0; i < downPassSeq.size(); i++) {
        double cstate = ostree->getNodeFromDownpassSeq(i)->getBrlen();
        getNodeFromDownpassSeq(i)->setTraitValue(cstate);

    }


}


double Tree::getTraitMaxTip(void)
{

    double ctrait = root->getTraitValue();
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getLfDesc() == NULL && (*i)->getTraitValue() > ctrait )
            ctrait = (*i)->getTraitValue();
    }
    return ctrait;
}

double Tree::getTraitMinTip(void)
{
    double ctrait = root->getTraitValue();
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getLfDesc() == NULL && (*i)->getTraitValue() < ctrait )
            ctrait = (*i)->getTraitValue();
    }
    return ctrait;

}

void Tree::printNodeBranchRates(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        cout << (*i)->getMeanSpeciationRate() << "\t" << (*i)->getNodeLambda() << endl;

}

void Tree::printNodeTraitRates(void)
{
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++)
        cout << (*i) << "\t" << (*i)->getMeanBeta() << endl;

}



void Tree::setMeanBranchTraitRates(void)
{

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        computeMeanTraitRatesByNode((*i));
        //cout  << (*i) << "\t" << (*i)->getMeanBeta() << "\t" << (*i)->getTraitBranchHistory()->getAncestralNodeEvent() << endl;

    }

}


/*
    This function is replicated for speciation and extinction as part of
    class node - it seems more efficient to put it here with class tree.


 */

void Tree::computeMeanTraitRatesByNode(Node* x)
{

    TraitBranchHistory* bh = x->getTraitBranchHistory();

    if (x->getAnc() != NULL) {
        // Only compute mean branch rate if node is NOT the root



        double rate = 0.0;
        int n_events = bh->getNumberOfBranchEvents();

        if (n_events == 0) {

            double t1 = x->getAnc()->getTime();
            double t2 = x->getTime();
            // times must be relative to event occurrence time:
            t1 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            t2 -= bh->getAncestralNodeEvent()->getAbsoluteTime();

            double zpar = bh->getAncestralNodeEvent()->getBetaShift();
            double beta0 = bh->getAncestralNodeEvent()->getBetaInit();

            if (zpar == 0)
                rate = beta0;
            else {
                rate = (beta0 / zpar) * ( exp(zpar * t2) - exp(zpar * t1));
                rate /= x->getBrlen();
            }

        } else {

            double tcheck = 0.0;
            double t1 = x->getAnc()->getTime();
            double t2 = bh->getEventByIndexPosition(0)->getAbsoluteTime();

            tcheck += (t2 - t1);

            // times must be relative to initial time of event
            t1 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            t2 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            double zpar = bh->getAncestralNodeEvent()->getBetaShift();
            double beta0 = bh->getAncestralNodeEvent()->getBetaInit();

            if (zpar == 0)
                rate = beta0 * (t2 - t1);
            else {
                //rate = frac * ((lam0/zpar) * ( exp(zpar*t2) - exp( zpar * t1)) / (t2-t1));
                // This can be simplified. All we have to do is sum the integrals over different
                //  rate-portions of the branch, then divide the result by the branch length

                rate = ((beta0 / zpar) * ( exp(zpar * t2) - exp( zpar * t1)));
            }

            //cout << x << "\t" << "1st event on branch: " << endl;
            //cout << "T1: " << t1 <<  "\tT2: " << t2 << endl;
            //cout << "rate pars: b0\t" << beta0 << "\tzpar:\t" << zpar << endl;



            for (int k = 1; k < n_events; k++) {

                //t1 = t2;
                t1 = 0.0;
                t2 = bh->getEventByIndexPosition(k)->getAbsoluteTime() -
                     bh->getEventByIndexPosition((k - 1))->getAbsoluteTime();
                zpar = bh->getEventByIndexPosition((k - 1))->getBetaShift();
                beta0 = bh->getEventByIndexPosition((k - 1))->getBetaInit();

                if (zpar == 0)
                    rate += beta0 * (t2 - t1);
                else
                    rate += (beta0 / zpar) * ( exp(zpar * t2) - exp(zpar * t1));
                //cout << k-1 <<  "\tt1: " <<  t1 << "\tt2: " << t2 <<  "\t" << t2 - t1 << "\t" << tcheck<<endl;

                tcheck += (t2 - t1);
            }

            //t1 = t2;
            t1 = 0.0;
            t2 = x->getTime() - bh->getEventByIndexPosition((n_events -
                    1))->getAbsoluteTime();

            zpar = bh->getNodeEvent()->getBetaShift();
            beta0 = bh->getNodeEvent()->getBetaInit();
            if (zpar == 0)
                rate += beta0 * (t2 - t1);
            else
                rate += (beta0 / zpar) * ( exp(zpar * t2) - exp(zpar * t1));
            tcheck += (t2 - t1);

            //cout << x << "\t2nd event on branch: " << endl;
            //cout << "T1: " << t1 <<  "\tT2: " << t2 << endl;
            //cout << "rate pars: b0\t" << beta0 << "\tzpar:\t" << zpar << endl;
            //cout << "Branch length: " << x->getBrlen() << endl << endl;

            // The overall mean rate across the branch:
            rate /= (x->getBrlen());

            //cout << "Rate: " << rate << endl;

        }
        x->setMeanBeta(rate);

    } else {
        // Node is root
        x->setMeanBeta((double)0.0);

    }

    // compute speciation rate at the focal node:
    double reltime = x->getTime() - bh->getNodeEvent()->getAbsoluteTime();
    double curBeta = bh->getNodeEvent()->getBetaInit() * exp((
                         reltime * bh->getNodeEvent()->getBetaShift()));


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
void Tree::initializeTraitValues(void)
{

    cout << "Setting initial trait values at internal nodes" << endl;

    // get min & max values:
    double mn = 0;
    double mx = 0;
    bool set = false;

    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {

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

    /*

     for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++){
     cout << (*i) << "\tFixed" << (*i)->getIsTraitFixed() << "\tValue: " << (*i)->getTraitValue() << endl;
     }
     */
    /*

     for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++){
     if ((*i)->getIsTraitFixed() == false){
     (*i)->setTraitValue( ranPtr->uniformRv(mn, mx) );
     }
     }
     */

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

        } else if (s2 > s1) {
            //x->setTraitValue(ranPtr->uniformRv(s2, s1));
            x->setTraitValue((s1 + s2) / (double)2);
        } else
            x->setTraitValue(s1);
    } else if (x->getIsTraitFixed() == false)
        x->setTraitValue(ranPtr->uniformRv(mn, mx));
    else {
        // Trait is fixed. Nothing to do.
    }
}


// Returns pointer to node of mrca of 2 taxa, with names
//  A and B.

Node* Tree::getNodeMRCA(string A, string B)
{
    //cout << "MRCA of " << A << "\t" << B << endl;

    Node* nodeA;
    Node* nodeB;
    bool Agood = false;
    bool Bgood = false;

    clearTempNodeArray();
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
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
        cout << "invalid nodes sent to Tree::getNodeMRCA(...)" << endl;
        cout << "\nEXITING WITH ERROR\n" << endl;
        exit(0);
    }

    passUpFillTempNodeArray(nodeA);
    //cout << _tempNodeSet.size() << endl;

    bool isFoundMRCA = false;
    while (!isFoundMRCA) {
        nodeB = nodeB->getAnc();
        //cout << nodeB << endl;
        if (_tempNodeSet.count(nodeB) > 0) {
            //cout << nodeB << " in common" << endl;
            break;
        }


    }
    //cout << endl << endl;
    //cout << nodeB << endl;
    //for (std::set<Node*>::iterator i = _tempNodeSet.begin(); i != _tempNodeSet.end(); i++)
    //  cout << "From A: " << (*i) << endl;
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

Node* Tree::getNodeByName(string A)
{

    Node* x = root;
    int count = 0;
    for (std::set<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        if ((*i)->getName() == A) {
            count++;
            x = (*i);
        }
    }
    if (count == 0) {
        cout << "Invalid node name: name not found in Tree:: getNodeByName" << endl;
        exit(0);
    } else if (count > 1) {
        cout << "Duplicate node names found in Tree:: getNodeByName" << endl;
        exit(0);
    } else {

    }

    return x;
}

















