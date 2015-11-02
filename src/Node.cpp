#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "Node.h"
#include "BranchHistory.h"
#include "SpExBranchEvent.h"


Node::Node()
{
    init();
}


Node::Node(int x)
{
    init(x);
}


void Node::init(int x)
{
    _lfDesc = NULL;
    _rtDesc = NULL;
    _anc = NULL;
    _name = "";
    _index = x;
    _time = 0.0;
    _brlen = 0.0;
    _isTip = false;
    _isExtant = false;
    _isConstant = false;
    _tipDescCount = 0;

    _mapStart = 0.0;
    _mapEnd = 0.0;

    _branchTime = 0.0;
    _cladeName = "";

    _history = new BranchHistory();

    // Compound Poisson stuff
    _meanSpeciationRate = 0.0;
    _meanExtinctionRate = 0.0;

    _nodeLambda = 0.0;
    _nodeMu = 0.0;

    // For phenotypes:
    //Phenotype * pheno = new Phenotype();
    _trait = 0.0;
    _meanBeta = 0;
    _isTraitFixed = 0;

    _ei = -1.0;
    _di = -1.0;
    _etip = -1.0;
    _nodeLikelihood = 0.0;

    _canHoldEvent = false;
    
    _eEnd = -1.0;
    _hasDownstreamRateShift = false;
    _inheritFromLeft = false;
}


// The equation for lambda through time, lam(t), is
//
//                 /  lam0 (2 - e^(-kt))    k > 0
//                |
//     lam(t) =  <    lam0 e^(kt)           k < 0
//                |
//                 \  lam0                  k = 0
//
// where lam0 is the initial lambda, k is the rate parameter, and t is time.
// The definite integral, named "Lam" here, is calculated as
//
//            t2
//           /
//     Lam = | lam(t) dt
//           /
//            t1
//
// where t1 and t2 are the times bounding the integral.
// Solving this definite integral results in
//
//             /  lam0 [2(t2 - t1) + (1/k)(e^(-k*t2) - e^(-k*t1))]   k > 0
//            |
//            |   lam0
//     Lam = <    ---- [e^(k*t2) - e^(k*t1)]                         k < 0
//            |    k
//            |
//             \  lam0 * (t2 - t1)                                   k = 0

double Node::integrateExponentialRateFunction
    (double init, double shift, double t1, double t2)
{
    if (shift < 0) {
        return (init / shift) * (std::exp(shift * t2) - std::exp(shift * t1));
    } else if (shift > 0) {
        return init * (2 * (t2 - t1) + (1.0 / shift) *
            (std::exp(-shift * t2) - std::exp(-shift * t1)));
    } else {
        return init * (t2 - t1);
    }
}


// The equation for lambda through time, lam(t), is
//
//                 /  lam0 (2 - e^(-kt))    k > 0
//                |
//     lam(t) =  <    lam0 e^(kt)           k < 0
//                |
//                 \  lam0                  k = 0

double Node::getExponentialRate(double init, double shift, double t)
{
    if (shift < 0) {
        return init * std::exp(shift * t);
    } else if (shift > 0) {
        return init * (2 - std::exp(-shift * t));
    } else {
        return init;
    }
}


int Node::getDescCount()
{
    int count = 0;
    if (getLfDesc() != NULL)
        count++;
    if (getRtDesc() != NULL)
        count++;
    return count;
}

Node* Node::getRandomLeftTipNode()
{


    if (getLfDesc() != NULL) {
        Node* xnode = getLfDesc();
        if (xnode->getLfDesc() == NULL)
            return xnode;
        else {
            while (xnode->getLfDesc() != NULL)
                xnode = xnode->getLfDesc();
            return xnode;
        }
    } else
        return NULL;
}

Node* Node::getRandomRightTipNode()
{

    if (getRtDesc() != NULL) {
        Node* xnode = getRtDesc();
        if (xnode->getLfDesc() == NULL)
            return xnode;
        else {
            while (xnode->getLfDesc() != NULL)
                xnode = xnode->getLfDesc();
            return xnode;
        }
    } else
        return NULL;
}


double Node::pathLengthToRoot()
{
    double length = getBrlen();

    Node* node = this;
    while (node->getAnc() != NULL) {
        node = node->getAnc();
        length += node->getBrlen();
    }

    return length;
}


/* MARCH 24 2012
 SET SPECIATION & EXTINCTION RATES BY NODE



 */

// If arg == void
//      COMPUTE rates
void Node::computeNodeBranchSpeciationParams()
{

    BranchHistory* bh = getBranchHistory();

    if (getAnc() != NULL) {
        SpExBranchEvent* ancestralEvent =
            static_cast<SpExBranchEvent*>(bh->getAncestralNodeEvent());

        double rate = 0.0;
        int n_events = bh->getNumberOfBranchEvents();

        if (n_events == 0) {

            double t1 = getAnc()->getTime();
            double t2 = getTime();

            // Times must be relative to event occurrence time
            t1 -= ancestralEvent->getAbsoluteTime();
            t2 -= ancestralEvent->getAbsoluteTime();

            double zpar = ancestralEvent->getLamShift();
            double lam0 = ancestralEvent->getLamInit();

            rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);
            rate /= getBrlen();

        } else {

            double tcheck = 0.0;
            double t1 = getAnc()->getTime();
            double t2 = bh->getEventByIndexPosition(0)->getAbsoluteTime();

            tcheck += (t2 - t1);

            // Times must be relative to initial time of event
            t1 -= ancestralEvent->getAbsoluteTime();
            t2 -= ancestralEvent->getAbsoluteTime();

            double zpar = ancestralEvent->getLamShift();
            double lam0 = ancestralEvent->getLamInit();

            rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);

            for (int k = 1; k < n_events; k++) {
                SpExBranchEvent* eventAtK = static_cast<SpExBranchEvent*>
                    (bh->getEventByIndexPosition(k));
                SpExBranchEvent* eventAtKMinus1 = static_cast<SpExBranchEvent*>
                    (bh->getEventByIndexPosition(k - 1));

                t1 = 0.0;
                t2 = eventAtK->getAbsoluteTime() -
                     eventAtKMinus1->getAbsoluteTime();
                zpar = eventAtKMinus1->getLamShift();
                lam0 = eventAtKMinus1->getLamInit();

                rate += integrateExponentialRateFunction(lam0, zpar, t1, t2);

                tcheck += (t2 - t1);
            }

            t1 = 0.0;
            t2 = getTime() -
                bh->getEventByIndexPosition((n_events - 1))->getAbsoluteTime();

            SpExBranchEvent* event = static_cast<SpExBranchEvent*>
                (bh->getNodeEvent());

            zpar = event->getLamShift();
            lam0 = event->getLamInit();

            rate += integrateExponentialRateFunction(lam0, zpar, t1, t2);

            tcheck += (t2 - t1);

            // The overall mean rate across the branch:
            rate /= (getBrlen());
        }

        setMeanSpeciationRate(rate); // extinction rate for branch set

    } else {
        // Node is root
        setMeanSpeciationRate(0.0);
    }

    // Compute speciation rate at the focal node
    SpExBranchEvent* event = static_cast<SpExBranchEvent*>(bh->getNodeEvent());
    double reltime = getTime() - event->getAbsoluteTime();
    double r_init = event->getLamInit();
    double r_shift = event->getLamShift();
    double curLam = getExponentialRate(r_init, r_shift, reltime);

#ifdef DEBUG_TIME_VARIABLE
    // Try setting node speciation rates equal to mean rate on
    // descendant branches, to see if the high-rate trap disappears
#else
    setNodeLambda(curLam); // speciation rate for node set
#endif
}



void Node::computeAndSetNodeSpeciationParams()
{
    BranchHistory* bh = getBranchHistory();
    // Compute speciation rate at the focal node
    SpExBranchEvent* event = static_cast<SpExBranchEvent*>(bh->getNodeEvent());
    double reltime = getTime() - event->getAbsoluteTime();
    double r_init = event->getLamInit();
    double r_shift = event->getLamShift();
    double curLam = getExponentialRate(r_init, r_shift, reltime);

    setNodeLambda(curLam); // speciation rate for node set

}

void Node::computeAndSetNodeExtinctionParams()
{

    BranchHistory* bh = getBranchHistory();
    // Compute extinction rate at the focal node:
    SpExBranchEvent* event = static_cast<SpExBranchEvent*>(bh->getNodeEvent());
    double reltime = getTime() - event->getAbsoluteTime();
    double r_init = event->getMuInit();
    double r_shift = event->getMuShift();
    double curMu = getExponentialRate(r_init, r_shift, reltime);

    setNodeMu(curMu); // extinction rate for node

}


void Node::computeNodeBranchExtinctionParams()
{

    BranchHistory* bh = getBranchHistory();

    if (getAnc() != NULL) {
        SpExBranchEvent* ancestralEvent =
            static_cast<SpExBranchEvent*>(bh->getAncestralNodeEvent());

        double rate = 0.0;
        int n_events = bh->getNumberOfBranchEvents();

        if (n_events == 0) {

            double t1 = getAnc()->getTime();
            double t2 = getTime();

            // Times must be relative to event occurrence time:
            t1 -= ancestralEvent->getAbsoluteTime();
            t2 -= ancestralEvent->getAbsoluteTime();

            double zpar = ancestralEvent->getMuShift();
            double r0 = ancestralEvent->getMuInit();

            rate = integrateExponentialRateFunction(r0, zpar, t1, t2);
            rate /= getBrlen();

        } else {

            double tcheck = 0.0;
            double t1 = getAnc()->getTime();
            double t2 = bh->getEventByIndexPosition(0)->getAbsoluteTime();

            tcheck += (t2 - t1);

            // Times must be relative to initial time of event
            t1 -= ancestralEvent->getAbsoluteTime();
            t2 -= ancestralEvent->getAbsoluteTime();
            double zpar = ancestralEvent->getMuShift();
            double r0 = ancestralEvent->getMuInit();

            rate = integrateExponentialRateFunction(r0, zpar, t1, t2);

            for (int k = 1; k < n_events; k++) {
                SpExBranchEvent* eventAtK = static_cast<SpExBranchEvent*>
                    (bh->getEventByIndexPosition(k));
                SpExBranchEvent* eventAtKMinus1 = static_cast<SpExBranchEvent*>
                    (bh->getEventByIndexPosition(k - 1));

                t1 = 0.0;
                t2 = eventAtK->getAbsoluteTime() -
                     eventAtKMinus1->getAbsoluteTime();
                zpar = eventAtKMinus1->getMuShift();
                r0 = eventAtKMinus1->getMuInit();

                rate += integrateExponentialRateFunction(r0, zpar, t1, t2);

                tcheck += (t2 - t1);
            }

            t1 = 0.0;
            t2 = getTime() -
                bh->getEventByIndexPosition((n_events - 1))->getAbsoluteTime();

            SpExBranchEvent* event = static_cast<SpExBranchEvent*>
                (bh->getNodeEvent());

            zpar = event->getMuShift();
            r0 = event->getMuInit();

            rate += integrateExponentialRateFunction(r0, zpar, t1, t2);

            tcheck += (t2 - t1);

            // The overall mean rate across the branch:
            rate /= (getBrlen());
        }

        setMeanExtinctionRate(rate); // extinction rate for branch set

    } else {
        // Node is root
        setMeanExtinctionRate(0.0);
    }

    // Compute extinction rate at the focal node:
    SpExBranchEvent* event = static_cast<SpExBranchEvent*>(bh->getNodeEvent());
    double reltime = getTime() - event->getAbsoluteTime();
    double r_init = event->getMuInit();
    double r_shift = event->getMuShift();
    double curMu = getExponentialRate(r_init, r_shift, reltime);

    setNodeMu(curMu); // extinction rate for node set
}


/*
 branchtime goes from 0 to brlen
 starting with t = 0 at ancestor


 */
double Node::getPointExtinction(double branchtime)
{
    BranchHistory* bh = getBranchHistory();
    double abstime = getTime() + getBrlen() - branchtime;
    double reltime = 0.0;
    double curMu = 0.0;
    if (bh->getNumberOfBranchEvents() == 0) {
        SpExBranchEvent* event =
            static_cast<SpExBranchEvent*>(bh->getNodeEvent());
        reltime = abstime - event->getAbsoluteTime();

        double shift = event->getMuShift();
        double init = event->getMuInit();

        curMu = getExponentialRate(init, shift, reltime);
    } else {
        // Multi-event scenario
        SpExBranchEvent* lastEvent =
            static_cast<SpExBranchEvent*>(bh->getAncestralNodeEvent());

        for (int i = 0; i < bh->getNumberOfBranchEvents(); i++) {
            SpExBranchEvent* temp_be =
                static_cast<SpExBranchEvent*>(bh->getEventByIndexPosition(i));
            if (temp_be->getAbsoluteTime() > abstime ) {
                break;
            }
            lastEvent = static_cast<SpExBranchEvent*>
                (bh->getEventByIndexPosition(i));
        }
        reltime = abstime - lastEvent->getAbsoluteTime();
        if (reltime < 0) {
            log(Error) << "Invalid time in Node::getPointExtinction().\n";
            std::exit(1);
        }

        double shift = lastEvent->getMuShift();
        double init = lastEvent->getMuInit();

        curMu = getExponentialRate(init, shift, reltime);
    }

    return curMu;
}

// These should be relative times.
double Node::computeSpeciationRateIntervalRelativeTime(double tstart,
        double tstop)
{
    // For FOSSIL process, do not check if tstop > getBrlen()
    //if ((tstart >= tstop) | (tstart < 0) | (tstop > getBrlen()) ) {
    
    if ((tstart >= tstop) | (tstart < 0) ) {
        std::cout << tstart << "\t" << tstop << std::endl;
        std::cout << "Invalid arguments to Node::computeSPeciationRateIntervalRelativeTime"
             << std::endl;
        throw;
    }
    BranchHistory* bh = getBranchHistory();

    double rate = 0.0;

    // COnvert start and stop times to absolute times...
    tstart += getAnc()->getTime();
    tstop += getAnc()->getTime();

    // This is required by the BranchHistory functions
    int n_events = bh->getNumberOfEventsOnInterval(tstart, tstop);
    //std::cout << "anctime: " << getAnc()->getTime() << "\tstart : " << tstart << "\tstop: " << tstop << std::endl;

    if (n_events == 0) {

        double t1 = tstart;
        double t2 = tstop;

        SpExBranchEvent* lastEvent =
            static_cast<SpExBranchEvent*>(bh->getLastEvent(tstart));

        t1 -= lastEvent->getAbsoluteTime();
        t2 -= lastEvent->getAbsoluteTime();

        double zpar = lastEvent->getLamShift();
        double lam0 = lastEvent->getLamInit();

        rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);
        rate /= (t2 - t1);

        if ((t1 < 0 ) | (t2 < t1)) {
            log(Error) << "Times are bad in "
                "Node::computeSpeciationRateIntervalRelTime.\n";
            std::exit(1);
        }

    } else {
        double tcheck = 0.0;
        double tabs1 = tstart;

        // Get next event from t1
        SpExBranchEvent* be =
            static_cast<SpExBranchEvent*>(bh->getLastEvent(tabs1));
        double tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();

        tcheck += (tabs2 - tabs1);

        // Times must be relative to initial time of event
        double trel1 = tstart - be->getAbsoluteTime();
        double trel2 = tabs2 - be->getAbsoluteTime();
        double zpar = be->getLamShift();
        double lam0 = be->getLamInit();

        rate = integrateExponentialRateFunction(lam0, zpar, trel1, trel2);

        be = static_cast<SpExBranchEvent*>(bh->getNextEvent(tabs1));

        for (int k = 1; k < n_events; k++) {

            tabs1 = be->getAbsoluteTime();
            tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();
            trel1 = 0.0;
            trel2 = tabs2 - tabs1;
            zpar = be->getLamShift();
            lam0 = be->getLamInit();

            rate += integrateExponentialRateFunction(lam0, zpar, trel1, trel2);

            tcheck += (trel2 - trel1);
            be = static_cast<SpExBranchEvent*>(bh->getNextEvent(tabs1));
        }

        trel1 = 0.0;
        trel2 = tstop - be->getAbsoluteTime();

        zpar = be->getLamShift();
        lam0 = be->getLamInit();

        rate += integrateExponentialRateFunction(lam0, zpar, trel1, trel2);

        tcheck += (trel2 - trel1);

        // The overall mean rate across the branch:
        rate /= tcheck;
    }

    return rate;
}

double Node::computeSpeciationRateIntervalRelativeTime(double t_init, double tstart,
                                                       double tstop)
{
    // For FOSSIL process, do not check if tstop > getBrlen()
    //if ((tstart >= tstop) | (tstart < 0) | (tstop > getBrlen()) ) {
    
    if ((tstart >= tstop) | (tstart < 0) ) {
        std::cout << "Invalid arguments to Node::computeSPeciationRateIntervalRelativeTime"
        << std::endl;
        throw;
    }
    BranchHistory* bh = getBranchHistory();
    
    double rate = 0.0;
    
    // COnvert start and stop times to absolute times...
    t_init += getAnc()->getTime();
    tstart += getAnc()->getTime();
    tstop += getAnc()->getTime();
        
    double t1 = tstart;
    double t2 = tstop;
        
    SpExBranchEvent* lastEvent = static_cast<SpExBranchEvent*>(bh->getLastEvent(t_init));
        
    t1 -= lastEvent->getAbsoluteTime();
    t2 -= lastEvent->getAbsoluteTime();
        
    double zpar = lastEvent->getLamShift();
    double lam0 = lastEvent->getLamInit();
        
    rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);
    rate /= (t2 - t1);
        
    if ((t1 < 0 ) | (t2 < t1)) {
        log(Error) << "Times are bad in "
        "Node::computeSpeciationRateIntervalRelTime.\n";
        std::exit(1);
    }
    
    return rate;
}



double Node::computeExtinctionRateIntervalRelativeTime(double t_init, double tstart, double tstop)
{
    
    if ((tstart >= tstop) | (tstart < 0) ) {
        std::cout << "Invalid arguments to Node::computeExtinctionRateIntervalRelativeTime"
        << std::endl;
        throw;
    }
    
    
    BranchHistory* bh = getBranchHistory();
    
    double rate = 0.0;
    
    t_init += getAnc()->getTime();
    tstart += getAnc()->getTime();
    tstop += getAnc()->getTime();
        
    double t1 = tstart;
    double t2 = tstop;
        
    SpExBranchEvent* lastEvent =
    static_cast<SpExBranchEvent*>(bh->getLastEvent(t_init));
        
    // Times must be relative to event occurrence time:
    t1 -= lastEvent->getAbsoluteTime();
    t2 -= lastEvent->getAbsoluteTime();
        
    double zpar = lastEvent->getMuShift();
    double mu0  = lastEvent->getMuInit();
        
    rate = integrateExponentialRateFunction(mu0, zpar, t1, t2);
    rate /= (t2 - t1);
        
    if ((t1 < 0 ) | (t2 < t1)) {
        log(Error) << "Times are bad in "
        "Node::computeExtinctionRateInterval.\n";
    std::exit(1);
    }
    
    return rate;
}



double Node::computeExtinctionRateIntervalRelativeTime(double tstart,
        double tstop)
{

    // For FOSSIL process, do not check if tstop > getBrlen()
    //if ((tstart >= tstop) | (tstart < 0) | (tstop > getBrlen()) ) {
    
    if ((tstart >= tstop) | (tstart < 0) ) {
        std::cout << "Invalid arguments to Node::computeExtinctionRateIntervalRelativeTime"
        << std::endl;
        throw;
    }


    BranchHistory* bh = getBranchHistory();

    double rate = 0.0;

    tstart += getAnc()->getTime();
    tstop += getAnc()->getTime();

    // This is required by the BranchHistory functions
    int n_events = bh->getNumberOfEventsOnInterval(tstart, tstop);

    if (n_events == 0) {

        double t1 = tstart;
        double t2 = tstop;

        SpExBranchEvent* lastEvent =
            static_cast<SpExBranchEvent*>(bh->getLastEvent(tstart));

        // Times must be relative to event occurrence time:
        t1 -= lastEvent->getAbsoluteTime();
        t2 -= lastEvent->getAbsoluteTime();

        double zpar = lastEvent->getMuShift();
        double mu0  = lastEvent->getMuInit();

        rate = integrateExponentialRateFunction(mu0, zpar, t1, t2);
        rate /= (t2 - t1);

        if ((t1 < 0 ) | (t2 < t1)) {
            log(Error) << "Times are bad in "
                "Node::computeExtinctionRateInterval.\n";
            std::exit(1);
        }
    } else {
        double tcheck = 0.0;
        double tabs1 = tstart;

        // Get next event from t1
        SpExBranchEvent* be =
            static_cast<SpExBranchEvent*>(bh->getLastEvent(tabs1));
        double tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();

        tcheck += (tabs2 - tabs1);

        // Times must be relative to initial time of event
        double trel1 = tstart - be->getAbsoluteTime();
        double trel2 = tabs2 - be->getAbsoluteTime();
        double zpar = be->getMuShift();
        double mu0 = be->getMuInit();

        rate = integrateExponentialRateFunction(mu0, zpar, trel1, trel2);

        be = static_cast<SpExBranchEvent*>(bh->getNextEvent(tabs1));
        for (int k = 1; k < n_events; k++) {
            tabs1 = be->getAbsoluteTime();
            tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();
            trel1 = 0.0;
            trel2 = tabs2 - tabs1;
            zpar = be->getMuShift();
            mu0 = be->getMuInit();

            rate += integrateExponentialRateFunction(mu0, zpar, trel1, trel2);

            tcheck += (trel2 - trel1);
            be = static_cast<SpExBranchEvent*>(bh->getNextEvent(tabs1));
        }

        trel1 = 0.0;
        trel2 = tstop - be->getAbsoluteTime();

        zpar = be->getMuShift();
        mu0 = be->getMuInit();

        rate += integrateExponentialRateFunction(mu0, zpar, trel1, trel2);

        tcheck += (trel2 - trel1);

        // The overall mean rate across the branch:
        rate /= tcheck;
    }

    return rate;
}


bool Node::isInternal()
{
    return (getLfDesc() != NULL) && (getRtDesc() != NULL);
}


std::string Node::getRandomRightDesc()
{
    Node* node = this;
    
    while (node->getRtDesc() != NULL){
        node = node->getRtDesc();
    }
    return node->getName();
}

std::string Node::getRandomLeftDesc()
{
    Node* node = this;
    
    while (node->getLfDesc() != NULL){
        node = node->getLfDesc();
    }
    return node->getName();

}





