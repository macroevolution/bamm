#include "TraitBranchEvent.h"

class Random;
class Node;
class Tree;


TraitBranchEvent::TraitBranchEvent(double beta, double shift,
    bool isTimeVariable, Node* x, Tree* tp, Random& random, double map) :
        BranchEvent(x, tp, random, map),
        _betaInit(beta), _betaShift(shift), _isTimeVariable(isTimeVariable)
{
}
