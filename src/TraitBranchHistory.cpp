
#include "TraitBranchHistory.h"



#include <iostream>
#include <stdio.h>

#include <stdlib.h>

using namespace std;

TraitBranchHistory::TraitBranchHistory(void)
{

    //cout << "BranchHistory ctor" << endl;
    nodeEvent = NULL;
    ancestralNodeEvent = NULL;


}


TraitBranchHistory::~TraitBranchHistory(void)
{



}

void TraitBranchHistory::printBranchHistory(void)
{

    int nEvents = eventsOnBranch.size();
    cout << "nodeEvent: " << nodeEvent << "\tancestorEvent: " << ancestralNodeEvent
         << endl;

    cout << "events on branch: " << nEvents << endl;
    if (nEvents > 0) {
        for (std::set<TraitBranchEvent*, comp_history>::iterator i =
                    eventsOnBranch.begin(); i != eventsOnBranch.end(); i++)
            cout << (*i) << "\t\t" << (*i)->getMapTime() << "\t" <<
                 (*i)->getAbsoluteTime() << endl;
    }

}


void TraitBranchHistory::reversePrintBranchHistory(void)
{

    cout << "Reverse printing events on branch: " << endl;

    std::set<TraitBranchEvent*>::iterator myIt = eventsOnBranch.end();
    myIt--;
    for ( ; myIt != eventsOnBranch.begin()--; myIt--)
        cout << "event at: " <<  (*myIt)->getMapTime() << endl;


}


TraitBranchEvent* TraitBranchHistory::getLastEvent(void)
{
    if (eventsOnBranch.size() == 0)
        return NULL;
    std::set<TraitBranchEvent*>::iterator myIt = eventsOnBranch.end();
    myIt--;
    return (*myIt);

}

// Return last event overall for a given event
//  assuming branch history is set correctly
TraitBranchEvent*  TraitBranchHistory::getLastEvent(TraitBranchEvent* x)
{

    TraitBranchEvent* theLastEvent = NULL;

    //int counter = 0;
    for (std::set<TraitBranchEvent*>::iterator myIt = eventsOnBranch.begin();
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
        cout << "problem in TraitBranchHistory::getLastEvent()" << endl;

    return   theLastEvent;

}


TraitBranchEvent* TraitBranchHistory::getEventByIndexPosition(int index)
{
    EventSetSizeType i = static_cast<EventSetSizeType>(index);
    if (i >= eventsOnBranch.size() ) {
        cout << "BranchHistory::getEventByIndexPosition error - accessing invalid event"
             << endl;
        return NULL;
    } else {

        std::set<TraitBranchEvent*>::iterator myIt = eventsOnBranch.begin();
        for (EventSetSizeType k = 0; k < i; k++)
            myIt++;
        return (*myIt);
    }

}






