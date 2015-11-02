#include <iostream>
#include <cstdlib>

#include "BranchEvent.h"
#include "Node.h"
#include "Tree.h"
#include "Random.h"
#include "Log.h"


BranchEvent::BranchEvent(Node* x, Tree* tp, Random& random, double map) :
    mapTime(map), nodeptr(x), treePtr(tp), _random(random),
    oldNodePtr(x), oldMapTime(map), _isEventTimeVariable(false)
{
    if (tp->getRoot() == x) {
        _absTime = 0.0;
    } else {
        _absTime = tp->getAbsoluteTimeFromMapTime(map);
    }
}


BranchEvent::~BranchEvent()
{
}


// If maptime is > than a.maptime, then event must have occurred earlier
// than a.maptime event; thus operator < would be true.
// For a->maptime > b->maptime, then a < b evaluates to true.

bool BranchEvent::operator<(const BranchEvent& a) const
{
    // TODO: This is counter-intuitive; explain
    if ( mapTime > a.mapTime ) {
        return true;
    } else {
        return false;
    }
}


// This (along with incrementMapPosition) implements a uniform proposal
// mechanism with reflecting boundaries at all major tree boundaries:
// at the root, the proposal is bounced back onto the parent branch.
// At the tips OR at invalid nodes (defined in tree constructor)

// This mechanism of local move involves actually moving the position of the
// event. If the proposed move is accepted, this is fine. If the proposed move
// is rejected, then the branchEvent will need to be reset to its original
// value. In practice, this should be done by:
//
//     oldMapPosition = getMapTime();
//     moveEventLocal();
//     calculate likelihood etc
//     if (accept)
//         event has already been updated.
//     if (reject)
//         setEventByMapPosition(oldMapPosition)

void BranchEvent::moveEventLocal(double stepsize)
{
    oldNodePtr = nodeptr;
    oldMapTime = mapTime;

    incrementMapPosition(stepsize);
}


void BranchEvent::incrementMapPosition(double ink)
{
    double mapstart = nodeptr->getMapStart();
    double mapend = nodeptr->getMapEnd();

    double temp = getMapTime() + ink;

    // Map end is ROOTWARDS
    // Map start is TIPWARDS

    if (temp > mapend) {
        // Goes to next most "ROOTWARDS" branch
        // BUT if there is no next branch (because ancestor of curr branch
        // is ROOT) go FORWARD down other descendant.

        ink = temp - mapend; // residual...

        if (getEventNode() == treePtr->getRoot()) {
            // Should never get here
            log(Error) << "Root problem in incrementMapPosition\n";
            std::exit(1);
        }

        if (getEventNode()->getAnc() == treePtr->getRoot()) {
            ink = -1 * ink;

            if ((treePtr->getRoot()->getLfDesc() == getEventNode()) &&
                    treePtr->getRoot()->getRtDesc()->getCanHoldEvent()) {
                setEventNode(treePtr->getRoot()->getRtDesc());

            } else if ((treePtr->getRoot()->getRtDesc() == getEventNode()) &&
                     treePtr->getRoot()->getLfDesc()->getCanHoldEvent()) {
                setEventNode(treePtr->getRoot()->getLfDesc());

            } else {
                // If you get here, either the right or left desc of the root
                // CANNOT hold event. Do nothing. Just go backwards and reflect
                // off of root.
            }

            setMapTime(getEventNode()->getMapEnd());
            incrementMapPosition(ink);

        } else {
            setEventNode(getEventNode()->getAnc());
            setMapTime(getEventNode()->getMapStart());
            incrementMapPosition(ink);
        }

    } else if (temp < mapstart) {

        ink = temp - mapstart; // ink is now negative...

        if (getEventNode()->getLfDesc() == NULL &&
                getEventNode()->getRtDesc() == NULL) {
            // At terminal branch:
            // Reflect
            ink  = -1 * ink; // make positive....
            setMapTime(mapstart);
            incrementMapPosition(ink);

        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == true) {
            // both desc branches valid
            double ran = _random.uniform();
            if (ran <= 0.5) {
                // left branch
                setEventNode(getEventNode()->getLfDesc());
                setMapTime(getEventNode()->getMapEnd());

            } else {
                setEventNode(getEventNode()->getRtDesc());
                setMapTime(getEventNode()->getMapEnd());
            }
            incrementMapPosition(ink);
        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == false &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == true) {
            setEventNode(getEventNode()->getRtDesc());
            setMapTime(getEventNode()->getMapEnd());
            incrementMapPosition(ink);

        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == true &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == false) {

            setEventNode(getEventNode()->getLfDesc());
            setMapTime(getEventNode()->getMapEnd());
            incrementMapPosition(ink);
        } else if (getEventNode()->getLfDesc()->getCanHoldEvent() == false &&
                   getEventNode()->getRtDesc()->getCanHoldEvent() == false) {
            // neither node can hold event: reflect
            ink  = -1 * ink; // make positive....

            // Reflect:
            setMapTime(mapstart);
            incrementMapPosition(ink);

        } else {
            log(Error) << "Problem in incrementMapPosition()\n";
            std::exit(1);
        }

    } else {
        setMapTime(temp);
    }

    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
}


void BranchEvent::moveEventGlobal()
{
    oldNodePtr = nodeptr;
    oldMapTime = mapTime;

    double aa = treePtr->getRoot()->getMapStart();
    double bb = treePtr->getTotalMapLength();
    double position = _random.uniform(aa, bb);
    setEventNode(treePtr->mapEventToTree(position));
    setMapTime(position);
    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
}


void BranchEvent::setEventByMapPosition(double x)
{
    setEventNode(treePtr->mapEventToTree(x));
    setMapTime(x);
    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
}


void BranchEvent::revertOldMapPosition()
{
    nodeptr = oldNodePtr;
    mapTime = oldMapTime;
    setAbsoluteTime(treePtr->getAbsoluteTimeFromMapTime(getMapTime()));
}
