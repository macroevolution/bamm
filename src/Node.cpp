#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "Node.h"
#include "BranchHistory.h"
#include "TraitBranchHistory.h"
//#include "MbRandom.h"

//#define _NEW_RATEFUNCTION


Node::Node(void)
{
    _lfDesc = NULL;
    _rtDesc = NULL;
    _anc = NULL;
    _name = "";
    _index = 0;
    _time = 0.0;
    _brlen = 0.0;
    //_descCount = NULL;
    _isTip = false;
    _isExtant = false;
    _isConstant = false;

    _mapStart = 0.0;
    _mapEnd = 0.0;

    _branchTime = 0.0;
    _cladeName = "";

    BranchHistory* bh = new BranchHistory();
    TraitBranchHistory* tbh = new TraitBranchHistory();

    history = bh;
    _traitHistory = tbh;

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

}

Node::Node(int x)
{
    _lfDesc = NULL;
    _rtDesc = NULL;
    _anc = NULL;
    _name = "";
    _index = x;
    _time = 0.0;
    _brlen = 0.0;
    //_descCount = NULL;
    _isExtant = false;
    _isTip = false;
    _isConstant = false;

    _mapStart = 0.0;
    _mapEnd = 0.0;

    _branchTime = 0.0;
    _cladeName = "";

    BranchHistory* bh = new BranchHistory();
    history = bh;

    TraitBranchHistory* tbh = new TraitBranchHistory();
    _traitHistory = tbh;

    _meanSpeciationRate = 0;
    _meanExtinctionRate = 0;

    _nodeLambda = 0.0;
    _nodeMu = 0.0;

    // For phenotypes:
    //Phenotype * pheno = new Phenotype();
    _trait = 0.0;
    _meanBeta = 0;
    _isTraitFixed = 0;

    //Speciation-extinction initial conditions
    _ei = -1.0;
    _di = -1.0;
    _etip = -1.0;
    _nodeLikelihood = 0.0;

    _canHoldEvent = false;

}


double Node::integrateExponentialRateFunction(double par_init, double shift, double t1, double t2)
{

    double x = 0.0;

#ifdef _NEW_RATEFUNCTION
    
    if (shift == 0){
        
        x += (t2 - t1) * par_init;
        
    }else if (shift > 0){
        
        shift *= -1.0;
        x += (2 * par_init * t2) - (2 * par_init * t1);
        x += (par_init / shift) * (exp(shift * t2) - exp(shift * t1));
        
    }else{
        
        x += (par_init / shift) * (exp(shift * t2) - exp(shift * t1));
    
    }

#else
    
    if (shift == 0){
        x += (t2 - t1) * par_init;   
    }else{
       x += (par_init / shift) * (exp(shift * t2) - exp(shift * t1));
    }
    
    
#endif
    
    return x;

}


double Node::getExponentialRate(double par_init, double shift, double tm)
{

    double rate = 0.0;
    
#ifdef _NEW_RATEFUNCTION
    if (shift > 0){
        shift *= -1.0;
        rate = 2*par_init - (par_init * exp(shift * tm));
    }else{
        rate = par_init * exp(shift * tm);
    }
    
#else
    rate = par_init * exp(shift * tm);

#endif
    
    return rate;

}


int Node::getDescCount(void)
{
    int count = 0;
    if (getLfDesc() != NULL)
        count++;
    if (getRtDesc() != NULL)
        count++;
    return count;
}

Node* Node::getRandomLeftTipNode(void)
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

Node* Node::getRandomRightTipNode(void)
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
void Node::computeNodeBranchSpeciationParams(void)
{

    BranchHistory* bh = getBranchHistory();

    if (getAnc() != NULL) {
        // Only compute mean branch rate if node is NOT the root



        double rate = 0.0;
        int n_events = bh->getNumberOfBranchEvents();

        if (n_events == 0) {

            double t1 = getAnc()->getTime();
            double t2 = getTime();
            // times must be relative to event occurrence time:
            t1 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            t2 -= bh->getAncestralNodeEvent()->getAbsoluteTime();

            double zpar = bh->getAncestralNodeEvent()->getLamShift();
            double lam0 = bh->getAncestralNodeEvent()->getLamInit();

            rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);
            rate /= getBrlen();


        } else {

            double tcheck = 0.0;
            double t1 = getAnc()->getTime();
            double t2 = bh->getEventByIndexPosition(0)->getAbsoluteTime();

            tcheck += (t2 - t1);

            // times must be relative to initial time of event
            t1 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            t2 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            double zpar = bh->getAncestralNodeEvent()->getLamShift();
            double lam0 = bh->getAncestralNodeEvent()->getLamInit();

            rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);

            for (int k = 1; k < n_events; k++) {

                //t1 = t2;
                t1 = 0.0;
                t2 = bh->getEventByIndexPosition(k)->getAbsoluteTime() -
                     bh->getEventByIndexPosition((k - 1))->getAbsoluteTime();
                zpar = bh->getEventByIndexPosition((k - 1))->getLamShift();
                lam0 = bh->getEventByIndexPosition((k - 1))->getLamInit();

                rate += integrateExponentialRateFunction(lam0, zpar, t1, t2);

                tcheck += (t2 - t1);
            }

            //t1 = t2;
            t1 = 0.0;
            t2 = getTime() - bh->getEventByIndexPosition((n_events - 1))->getAbsoluteTime();

            zpar = bh->getNodeEvent()->getLamShift();
            lam0 = bh->getNodeEvent()->getLamInit();
            
            rate += integrateExponentialRateFunction(lam0, zpar, t1, t2);

            tcheck += (t2 - t1);

            // The overall mean rate across the branch:
            rate /= (getBrlen());
        }

        setMeanSpeciationRate(rate); // extinction rate for branch set


    } else {
        // Node is root
        setMeanSpeciationRate((double)0.0);

    }

    // compute speciation rate at the focal node:
    double reltime = getTime() - bh->getNodeEvent()->getAbsoluteTime();
    double r_init = bh->getNodeEvent()->getLamInit();
    double r_shift = bh->getNodeEvent()->getLamShift();
    
    
    double curLam = getExponentialRate(r_init, r_shift , reltime);

#ifdef DEBUG_TIME_VARIABLE

    // Try setting node speciation rates equal to mean rate on descendant branches, to see if the
    //  high-rate trap disappears.



#else

    setNodeLambda(curLam); // speciation rate for node set

#endif


}




void Node::computeNodeBranchExtinctionParams(void)
{

    BranchHistory* bh = getBranchHistory();

    if (getAnc() != NULL) {
        // Only compute mean branch rate if node is NOT the root



        double rate = 0.0;
        int n_events = bh->getNumberOfBranchEvents();

        if (n_events == 0) {

            double t1 = getAnc()->getTime();
            double t2 = getTime();
            // times must be relative to event occurrence time:
            t1 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            t2 -= bh->getAncestralNodeEvent()->getAbsoluteTime();

            double zpar = bh->getAncestralNodeEvent()->getMuShift();
            double r0 = bh->getAncestralNodeEvent()->getMuInit();

            
            rate = integrateExponentialRateFunction(r0, zpar, t1, t2);
            rate /= getBrlen();


        } else {

            double tcheck = 0.0;
            double t1 = getAnc()->getTime();
            double t2 = bh->getEventByIndexPosition(0)->getAbsoluteTime();
            //std::cout << "Premod: t1: " << t1 << "\tt2: " << t2 << "\t" << t2 - t1 << std::endl;

            tcheck += (t2 - t1);

            // times must be relative to initial time of event
            t1 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            t2 -= bh->getAncestralNodeEvent()->getAbsoluteTime();
            double zpar = bh->getAncestralNodeEvent()->getMuShift();
            double r0 = bh->getAncestralNodeEvent()->getMuInit();

            rate = integrateExponentialRateFunction(r0, zpar, t1, t2);

            for (int k = 1; k < n_events; k++) {

                //t1 = t2;
                t1 = 0.0;
                t2 = bh->getEventByIndexPosition(k)->getAbsoluteTime() -
                     bh->getEventByIndexPosition((k - 1))->getAbsoluteTime();
                zpar = bh->getEventByIndexPosition((k - 1))->getMuShift();
                r0 = bh->getEventByIndexPosition((k - 1))->getMuInit();

                rate += integrateExponentialRateFunction(r0, zpar, t1, t2);

                tcheck += (t2 - t1);
            
            }

            //t1 = t2;
            t1 = 0.0;
            t2 = getTime() - bh->getEventByIndexPosition((n_events - 1))->getAbsoluteTime();

            zpar = bh->getNodeEvent()->getMuShift();
            r0 = bh->getNodeEvent()->getMuInit();
            
            rate += integrateExponentialRateFunction(r0, zpar, t1, t2);
    
            tcheck += (t2 - t1);
    
            // The overall mean rate across the branch:
            rate /= (getBrlen());
        }

        setMeanExtinctionRate(rate); // extinction rate for branch set

    } else {
        // Node is root
        setMeanExtinctionRate((double)0.0);

    }

    // compute extinction rate at the focal node:
    double reltime = getTime() - bh->getNodeEvent()->getAbsoluteTime();
    
    double r_init = bh->getNodeEvent()->getMuInit();
    double r_shift = bh->getNodeEvent()->getMuShift();
    double curMu = getExponentialRate(r_init, r_shift , reltime);


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
        reltime = abstime - bh->getNodeEvent()->getAbsoluteTime();
        
        double shift = bh->getNodeEvent()->getMuShift();
        double init = bh->getNodeEvent()->getMuInit();
        
        curMu = getExponentialRate(init, shift, reltime);

    } else {
        //multi-event scenario
        BranchEvent* lastEvent = bh->getAncestralNodeEvent();
        for (int i = 0; i < bh->getNumberOfBranchEvents(); i++) {
            BranchEvent* temp_be = bh->getEventByIndexPosition(i);
            if (temp_be->getAbsoluteTime() > abstime )
                break;
            lastEvent = bh->getEventByIndexPosition(i);
        }
        reltime = abstime - lastEvent->getAbsoluteTime();
        if (reltime < 0) {
            std::cout << "Invalid time in Node::getPointExtinction() " << std::endl;
            throw;
        }
        
        double shift = lastEvent->getMuShift();
        double init = lastEvent->getMuInit();
        
        curMu = getExponentialRate(init, shift, reltime);
        
    }


    return curMu;
}



double Node::computeSpeciationRateIntervalRelativeTime(double tstart,
        double tstop)
{

    if ((tstart >= tstop) | (tstart < 0) | (tstop > getBrlen()) ) {
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
        
        t1 -= bh->getLastEvent(tstart)->getAbsoluteTime();
        t2 -= bh->getLastEvent(tstart)->getAbsoluteTime();

        double zpar = bh->getLastEvent(tstart)->getLamShift();
        double lam0 = bh->getLastEvent(tstart)->getLamInit();

        rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);
        rate /= (t2 - t1);
        
        if ((t1 < 0 ) | (t2 < t1)) {
            std::cout << "error in Node::computeSpeciationRateIntervalRelTime - times are bad...\n" <<
                 std::endl;
            throw;
        }

    } else {
        double tcheck = 0.0;
        double tabs1 = tstart;
        // get next event from t1:
        BranchEvent* be = bh->getLastEvent(tabs1);
        double tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();

        tcheck += (tabs2 - tabs1);

        // times must be relative to initial time of event
        double trel1 = tstart - be->getAbsoluteTime();
        double trel2 = tabs2 - be->getAbsoluteTime();
        double zpar = be->getLamShift();
        double lam0 = be->getLamInit();

        rate = integrateExponentialRateFunction(lam0, zpar, trel1, trel2);
        
        be = bh->getNextEvent(tabs1);
        
        for (int k = 1; k < n_events; k++) {

            tabs1 = be->getAbsoluteTime();
            tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();
            trel1 = 0.0;
            trel2 = tabs2 - tabs1;
            zpar = be->getLamShift();
            lam0 = be->getLamInit();
            
            rate += integrateExponentialRateFunction(lam0, zpar, trel1, trel2);

            tcheck += (trel2 - trel1);
            be = bh->getNextEvent(tabs1);
        }

        //t1 = t2;
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


// What is this function doing? Where is it called?

double Node::computeSpeciationRateIntervalAbsoluteTime(double tstart,
        double tstop)
{

    if ((tstart >= tstop) | (tstart < getAnc()->getTime()) | (tstop > getTime()) ) {
        std::cout << "Invalid arguments to Node::computeSPeciationRateIntervalAbsoluteTime"
             << std::endl;
        throw;
    }

    BranchHistory* bh = getBranchHistory();

    double rate = 0.0;
    //std::cout << "at start: tstart: " << tstart << "\ttstop: " << tstop << std::endl;

    // COnvert start and stop times to absolute times...
    tstart += getAnc()->getTime();
    tstop += getAnc()->getTime();

    // This is required by the BranchHistory functions
    int n_events = bh->getNumberOfEventsOnInterval(tstart, tstop);
    //std::cout << "anctime: " << getAnc()->getTime() << "\tstart : " << tstart << "\tstop: " << tstop << std::endl;

    if (n_events == 0) {

        double t1 = tstart;
        double t2 = tstop;
        // times must be relative to event occurrence time:
        //std::cout << "LE time: " << bh->getLastEvent(tstart)->getAbsoluteTime() << std::endl;

        t1 -= bh->getLastEvent(tstart)->getAbsoluteTime();
        t2 -= bh->getLastEvent(tstart)->getAbsoluteTime();

        double zpar = bh->getLastEvent(tstart)->getLamShift();
        double lam0 = bh->getLastEvent(tstart)->getLamInit();

        //std::cout << "z: " << zpar << "\tlam0: " << lam0 << "\tt1: " << t1 << "\tt2: " << t2 << std::endl;

        rate = integrateExponentialRateFunction(lam0, zpar, t1, t2);
        rate /= (t2 - t1);

        //std::cout << "Rate: " <<  rate << std::endl;
        if ((t1 < 0 ) | (t2 < t1)) {
            std::cout << "error in Node::computeSpeciationRateIntervalAbsTime - times are bad...\n" <<
                 std::endl;
            throw;
        }

    } else {
        double tcheck = 0.0;
        double tabs1 = tstart;
        // get next event from t1:
        BranchEvent* be = bh->getLastEvent(tabs1);
        double tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();

        tcheck += (tabs2 - tabs1);

        // times must be relative to initial time of event
        double trel1 = tstart - be->getAbsoluteTime();
        double trel2 = tabs2 - be->getAbsoluteTime();
        double zpar = be->getLamShift();
        double lam0 = be->getLamInit();

        
        rate = integrateExponentialRateFunction(lam0, zpar, trel1, trel2);
    
        be = bh->getNextEvent(tabs1);
        for (int k = 1; k < n_events; k++) {

            tabs1 = be->getAbsoluteTime();
            tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();
            trel1 = 0.0;
            trel2 = tabs2 - tabs1;
            zpar = be->getLamShift();
            lam0 = be->getLamInit();
            
            rate += integrateExponentialRateFunction(lam0, zpar, trel1, trel2);

            tcheck += (trel2 - trel1);
            be = bh->getNextEvent(tabs1);
        }

        //t1 = t2;
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

double Node::computeExtinctionRateIntervalRelativeTime(double tstart,
        double tstop)
{

    if ((tstart >= tstop) | (tstart < 0) | (tstop > getBrlen()) ) {
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
        // times must be relative to event occurrence time:

        t1 -= bh->getLastEvent(tstart)->getAbsoluteTime();
        t2 -= bh->getLastEvent(tstart)->getAbsoluteTime();

        double zpar = bh->getLastEvent(tstart)->getMuShift();
        double mu0  = bh->getLastEvent(tstart)->getMuInit();

        rate = integrateExponentialRateFunction(mu0, zpar, t1, t2);
        rate /= (t2 - t1);

        if ((t1 < 0 ) | (t2 < t1)) {
            std::cout << "error in Node::computeExtinctionRateInterval - times are bad...\n" <<
                 std::endl;
            throw;
        }

    } else {
        double tcheck = 0.0;
        double tabs1 = tstart;
        // get next event from t1:
        BranchEvent* be = bh->getLastEvent(tabs1);
        double tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();

        tcheck += (tabs2 - tabs1);

        // times must be relative to initial time of event
        double trel1 = tstart - be->getAbsoluteTime();
        double trel2 = tabs2 - be->getAbsoluteTime();
        double zpar = be->getMuShift();
        double mu0 = be->getMuInit();

        rate = integrateExponentialRateFunction(mu0, zpar, trel1, trel2);
    
        be = bh->getNextEvent(tabs1);
        for (int k = 1; k < n_events; k++) {

            tabs1 = be->getAbsoluteTime();
            tabs2 = bh->getNextEvent(tabs1)->getAbsoluteTime();
            trel1 = 0.0;
            trel2 = tabs2 - tabs1;
            zpar = be->getMuShift();
            mu0 = be->getMuInit();
            
            rate += integrateExponentialRateFunction(mu0, zpar, trel1, trel2);

            tcheck += (trel2 - trel1);
            be = bh->getNextEvent(tabs1);
        }

        //t1 = t2;
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
