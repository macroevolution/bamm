#include "SpExBranchEvent.h"

class Node;
class Tree;
class MbRandom;


SpExBranchEvent::SpExBranchEvent(double speciation, double lamshift,
        double extinction, double mushift,
        Node* x, Tree* tp, MbRandom* rp, double map) :
    BranchEvent(x, tp, rp, map), _lamInit(speciation), _lamShift(lamshift),
        _muInit(extinction), _muShift(mushift)
{
}
