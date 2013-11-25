#include "TraitBranchHistory.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>


TraitBranchHistory::TraitBranchHistory(void)
{

    //std::cout << "BranchHistory ctor" << std::endl;
    nodeEvent = NULL;
    ancestralNodeEvent = NULL;


}


TraitBranchHistory::~TraitBranchHistory(void)
{



}

void TraitBranchHistory::printBranchHistory(void)
{

    int nEvents = (int)eventsOnBranch.size();
    std::cout << "nodeEvent: " << nodeEvent << "\tancestorEvent: " << ancestralNodeEvent
         << std::endl;

    std::cout << "events on branch: " << nEvents << std::endl;
    if (nEvents > 0) {
        for (EventSet::iterator i = eventsOnBranch.begin();
                i != eventsOnBranch.end(); i++)
            std::cout << (*i) << "\t\t" << (*i)->getMapTime() << "\t" <<
                 (*i)->getAbsoluteTime() << std::endl;
    }

}


void TraitBranchHistory::reversePrintBranchHistory(void)
{

    std::cout << "Reverse printing events on branch: " << std::endl;

    std::set<TraitBranchEvent*>::iterator myIt = eventsOnBranch.end();
    myIt--;
    for ( ; myIt != eventsOnBranch.begin()--; myIt--)
        std::cout << "event at: " <<  (*myIt)->getMapTime() << std::endl;


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
        //std::cout << "count: " << counter << std::endl;
    }

    if (theLastEvent == NULL)
        std::cout << "problem in TraitBranchHistory::getLastEvent()" << std::endl;

    return   theLastEvent;

}


TraitBranchEvent* TraitBranchHistory::getEventByIndexPosition(int index)
{
    EventSetSizeType i = static_cast<EventSetSizeType>(index);
    if (i >= eventsOnBranch.size() ) {
        std::cout << "BranchHistory::getEventByIndexPosition error - accessing invalid event"
             << std::endl;
        return NULL;
    } else {

        std::set<TraitBranchEvent*>::iterator myIt = eventsOnBranch.begin();
        for (EventSetSizeType k = 0; k < i; k++)
            myIt++;
        return (*myIt);
    }

}






