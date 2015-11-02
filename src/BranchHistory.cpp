#include "BranchHistory.h"
#include "BranchEvent.h"
#include "Log.h"

#include <cstdlib>


BranchHistory::BranchHistory() : _nodeEvent(NULL), _ancestralNodeEvent(NULL)
{
}


BranchEvent* BranchHistory::getLastEvent()
{
    if (_eventsOnBranch.size() == 0) {
        return NULL;
    }

    EventSet::reverse_iterator it = _eventsOnBranch.rbegin();
    return *it;
}


// Return the last event overall for a given event,
// assuming branch history is set correctly.
BranchEvent* BranchHistory::getLastEvent(BranchEvent* x)
{
    BranchEvent* theLastEvent = NULL;

    EventSet::iterator it;
    for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end(); ++it) {
        if (*it == x) {
            if (it == _eventsOnBranch.begin()) {
                theLastEvent = getAncestralNodeEvent();
            } else {
                theLastEvent = *(--it);
            }

            break;
        }
    }

    if (theLastEvent == NULL) {
        log(Warning) << "Problem in BranchHistory::getLastEvent()\n";
    }

    return theLastEvent;
}


// Returns most recent upstream event from a given absolute time,
// e.g., the most "rootward" event
BranchEvent* BranchHistory::getLastEvent(double ttime)
{
    BranchEvent* event = getAncestralNodeEvent();

    if (getNumberOfBranchEvents() > 0) {
        if (ttime > ((*_eventsOnBranch.begin())->getAbsoluteTime()) ) {
            if (_eventsOnBranch.size() == 1) {
                event = (*_eventsOnBranch.begin());
            } else {
                EventSet::iterator it;
                for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end();
                        ++it) {
                    if ((*it)->getAbsoluteTime() < ttime) {
                        event = *it;
                    } else {
                        break;
                    }
                }
            }
        } else {
            // Do nothing. Ancestral node event should be returned,
            // because ttime is BEFORE any events on branch
        }
    }

    return event;
}


// Returns most recent downstream event from a given absolute time,
// e.g.. the most "tipward" event.
// Does not check for range violation, e.g., if ttime is not in the focal branch
// If no events on branch, this will just be the EventNode.
BranchEvent* BranchHistory::getNextEvent(double ttime)
{
    BranchEvent* event = getNodeEvent();

    if (getNumberOfBranchEvents() > 0) {
        if (ttime < ((*(--_eventsOnBranch.end()))->getAbsoluteTime()) ) {
            if (_eventsOnBranch.size() == 1) {
                event = (*_eventsOnBranch.begin());
            } else {
                EventSet::iterator it;
                for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end();
                        ++it) {
                    if ((*it)->getAbsoluteTime() > ttime) {
                        event = *it;
                        break;
                    }
                }
            }
        } else {
            // Do nothing. EventNode should be returned
            // because ttime is after all events on branch
        }
    }

    return event;
}


// t1, t2 must be absolute time
int BranchHistory::getNumberOfEventsOnInterval(double t1, double t2)
{
    int n_events = 0;

    EventSet::iterator it;
    for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end(); ++it) {
        double atime = (*it)->getAbsoluteTime();
        //if ((atime > t1) && (atime <= t2)) {
        if ((atime >= t1) && (atime < t2)) {
            n_events++;
        }
    }

    return n_events;
}


BranchEvent* BranchHistory::getEventByIndexPosition(int index)
{
    EventSetSizeType i = static_cast<EventSetSizeType>(index);
    if (i < _eventsOnBranch.size() ) {
        EventSet::iterator it = _eventsOnBranch.begin();
        for (EventSetSizeType k = 0; k < i; k++)
            it++;
        return *it;
    } else {
        log(Error) << "BranchHistory::getEventByIndexPosition: "
                   << "accessing invalid event\n";
        std::exit(1);
    }

}


void BranchHistory::printBranchHistory()
{
    log() << "Node Event: " << _nodeEvent << "\t"
          << "Ancestor Event: " << _ancestralNodeEvent << "\n";

    EventSetSizeType numEvents = _eventsOnBranch.size();
    log() << "Number of events on branch: " << numEvents << "\n";

    EventSet::iterator it;
    for (it = _eventsOnBranch.begin(); it != _eventsOnBranch.end(); ++it) {
        printEvent(*it);
    }
}


void BranchHistory::reversePrintBranchHistory()
{
    EventSet::reverse_iterator it;
    for (it = _eventsOnBranch.rbegin(); it != _eventsOnBranch.rend(); ++it) {
        printEvent(*it);
    }
}


void BranchHistory::printEvent(BranchEvent* event)
{
    log() << event << "\t\t"
          << event->getMapTime() << "\t\t"
          << event->getAbsoluteTime() << "\n";
}


void BranchHistory::setNodeEvent(BranchEvent* x)
{
    _nodeEvent = x;
}


BranchEvent* BranchHistory::getNodeEvent()
{
    return _nodeEvent;
}


void BranchHistory::setAncestralNodeEvent(BranchEvent* x)
{
    _ancestralNodeEvent = x;
}


BranchEvent* BranchHistory::getAncestralNodeEvent()
{
    return _ancestralNodeEvent;
}


void BranchHistory::popEventOffBranchHistory(BranchEvent* x)
{
    _eventsOnBranch.erase(x);
}


void BranchHistory::addEventToBranchHistory(BranchEvent* x)
{
    _eventsOnBranch.insert(x);
}


int BranchHistory::getNumberOfBranchEvents()
{
    return (int)_eventsOnBranch.size();
}
