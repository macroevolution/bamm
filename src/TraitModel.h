#ifndef TRAIT_MODEL_H
#define TRAIT_MODEL_H

#include "Model.h"
#include "BetaInitProposal.h"
#include "BetaShiftProposal.h"
#include "NodeStateProposal.h"

#include <iosfwd>

class Node;
class MbRandom;
class Settings;
class Prior;
class BranchEvent;
class Proposal;


class TraitModel : public Model
{

public:

    TraitModel(MbRandom* rng, Settings* settings, Prior* prior);

    virtual double computeLogLikelihood();
    virtual double computeTriadLikelihoodTraits(Node* x);

    virtual double computeLogPrior();

private:

    virtual void readModelSpecificParameters(std::ifstream& inputFile);
    virtual void setRootEventWithReadParameters();

    void initializeSpecificUpdateWeights();
    virtual Proposal* getSpecificProposal(int parameter);

    virtual BranchEvent* newBranchEventWithReadParameters(Node* x, double time);
    virtual BranchEvent* newBranchEventWithRandomParameters(double x);
    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    virtual void setMeanBranchParameters();
    virtual void setDeletedEventParameters(BranchEvent* be);

    virtual double calculateLogQRatioJump();

    virtual void getSpecificEventDataString
        (std::stringstream& ss, BranchEvent* event);

    BetaInitProposal _betaInitProposal;
    BetaShiftProposal _betaShiftProposal;
    NodeStateProposal _nodeStateProposal;

    bool _sampleFromPriorOnly;

    double _lastDeletedEventBetaInit;;
    double _lastDeletedEventBetaShift;

    // Here are several variables that track the previous
    // state. At some point, these should have their own class

    // Ultimately, initializations should be handled in TraitModel
    // and Model classes
    // NOT in class TREE!!
    void   initializeTraitParamsForNodes();

    double _readBetaInit;
    double _readBetaShift;
};


#endif
