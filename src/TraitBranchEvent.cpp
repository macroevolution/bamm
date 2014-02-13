#include "TraitBranchEvent.h"

class MbRandom;
class Node;
class Tree;


TraitBranchEvent::TraitBranchEvent(double beta, double shift,
        Node* x, Tree* tp, MbRandom* rp, double map) :
    BranchEvent(x, tp, rp, map), _betaInit(beta), _betaShift(shift)
{
}
