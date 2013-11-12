/*
 *  BranchHistory.cpp
 *  rateshift
 *
 *  Created by Dan Rabosky on 12/14/11.
 *
 *
 */

#include <iostream>
#include <stdio.h>

#include <stdlib.h>

#include "BranchHistory.h"

using namespace std;

BranchHistory::BranchHistory(void)
{

    //cout << "BranchHistory ctor" << endl;
    nodeEvent = NULL;
    ancestralNodeEvent = NULL;


}


BranchHistory::~BranchHistory(void)
{



}

void BranchHistory::printBranchHistory(void)
{

    int nEvents = eventsOnBranch.size();
    cout << "nodeEvent: " << nodeEvent << "\tancestorEvent: " << ancestralNodeEvent
         << endl;

    cout << "events on branch: " << nEvents << endl;
    if (nEvents > 0) {
        for (std::set<BranchEvent*, comp_history>::iterator i = eventsOnBranch.begin();
                i != eventsOnBranch.end(); i++)
            cout << (*i) << "\t\t" << (*i)->getMapTime() << "\t" <<
                 (*i)->getAbsoluteTime() << endl;
    }

}


void BranchHistory::reversePrintBranchHistory(void)
{

    cout << "Reverse printing events on branch: " << endl;

    std::set<BranchEvent*>::iterator myIt = eventsOnBranch.end();
    myIt--;
    for ( ; myIt != eventsOnBranch.begin()--; myIt--)
        cout << "event at: " <<  (*myIt)->getMapTime() << endl;


}


BranchEvent* BranchHistory::getLastEvent(void)
{
    if (eventsOnBranch.size() == 0)
        return NULL;
    std::set<BranchEvent*>::iterator myIt = eventsOnBranch.end();
    myIt--;
    return (*myIt);

}

// Return last event overall for a given event
//  assuming branch history is set correctly
BranchEvent* BranchHistory::getLastEvent(BranchEvent* x)
{

    BranchEvent* theLastEvent = NULL;

    //int counter = 0;
    for (std::set<BranchEvent*>::iterator myIt = eventsOnBranch.begin();
            myIt != eventsOnBranch.end(); myIt++) {
        //counter++;
        if (((*myIt) == x) && myIt == eventsOnBranch.begin()) {

            theLastEvent = getAncestralNodeEvent();
            break;
        } else if ( (*myIt) == x ) {
            myIt--;
            theLastEvent = (*myIt);
            break;
        } else {

        }
        //cout << "count: " << counter << endl;
    }

    if (theLastEvent == NULL)
        cout << "problem in BranchHistory::getLastEvent()" << endl;

    return   theLastEvent;

}

// Returns most recent UPSTREAM event from a given absolute time
// E.g.. the most "rootward" event


BranchEvent* BranchHistory::getLastEvent(double ttime)
{


    BranchEvent* be = getAncestralNodeEvent();

    if (getNumberOfBranchEvents() > 0) {

        if (ttime > ((*eventsOnBranch.begin())->getAbsoluteTime()) ) {

            if (eventsOnBranch.size() == 1)
                be = (*eventsOnBranch.begin());
            else {
                for (std::set<BranchEvent*>::iterator i = eventsOnBranch.begin();
                        i != eventsOnBranch.end(); i++) {
                    if ((*i)->getAbsoluteTime() < ttime)
                        be = (*i);
                    else
                        break;
                }
            }
        } else {
            // do nothing. Ancestral node event should be returned,
            // because ttime is BEFORE any events on branch
        }


    }
    return be;
}

// Returns most recent DOWNSTREAM event from a given absolute time
// E.g.. the most "tipward" event
// DOES NOT check for range violation, eg if ttime is not in the focal branch
// If no events on branch, this will just be the EventNode.

BranchEvent* BranchHistory::getNextEvent(double ttime)
{


    BranchEvent* be = getNodeEvent();

    if (getNumberOfBranchEvents() > 0) {
        if (ttime < ((*(--eventsOnBranch.end()))->getAbsoluteTime()) ) {
            if (eventsOnBranch.size() == 1)
                be = (*eventsOnBranch.begin());
            else {
                for (std::set<BranchEvent*>::iterator i = eventsOnBranch.begin();
                        i != eventsOnBranch.end(); i++) {
                    if ((*i)->getAbsoluteTime() > ttime) {
                        be = (*i);
                        break;
                    } else {

                    }
                }
            }
        } else {
            // do nothing. EventNode should be returned
            // because ttime is AFTER all events on branch
        }


    }

    return be;
}

// t1, t2 must be ABSOLUTE time
int BranchHistory::getNumberOfEventsOnInterval(double t1, double t2)
{
    int n_events = 0;
    if (eventsOnBranch.size() > 0) {
        for (std::set<BranchEvent*>::iterator i = eventsOnBranch.begin();
                i != eventsOnBranch.end(); i++) {
            double atime = (*i)->getAbsoluteTime();
            if ((atime > t1) & (atime <= t2) )
                n_events++;
        }
    } else {
        // continue
    }
    return n_events;
}



BranchEvent* BranchHistory::getEventByIndexPosition(int index)
{
    EventSetSizeType i = static_cast<EventSetSizeType>(index);
    if (i >= eventsOnBranch.size() ) {
        cout << "BranchHistory::getEventByIndexPosition error - accessing invalid event"
             << endl;
        return NULL;
    } else {

        std::set<BranchEvent*>::iterator myIt = eventsOnBranch.begin();
        for (EventSetSizeType k = 0; k < i; k++)
            myIt++;
        return (*myIt);
    }

}






