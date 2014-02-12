#ifndef TRAIT_MODEL_H
#define TRAIT_MODEL_H

#include "Model.h"

#include <set>
#include <vector>
#include <sstream>

#include "BranchEvent.h"

// Forward declarations
class Tree;
class Node;
class MbRandom;
class Settings;
class Prior;


class TraitModel : public Model
{

public:

    TraitModel(MbRandom* rng, Tree* tree, Settings* settings, Prior* prior);
    virtual ~TraitModel();

    virtual double computeLogLikelihood();
    virtual double computeTriadLikelihoodTraits(Node* x);

    virtual double computeLogPrior();

    void updateBetaMH();
    void updateNodeStateMH();
    void updateNodeStateMH(Node* xnode);
    void updateBetaShiftMH();
    void setMinMaxTraitPriors();

    double      getLastLH();

private:

    virtual void readModelSpecificParameters(std::ifstream& inputFile);
    virtual void setRootEventWithReadParameters();
    virtual BranchEvent* newBranchEventWithReadParameters(Node* x, double time);
    virtual void setMeanBranchParameters();

    virtual BranchEvent* newBranchEventWithRandomParameters(double x);

    virtual void setDeletedEventParameters(BranchEvent* be);
    virtual double calculateLogQRatioJump();

    virtual BranchEvent* newBranchEventFromLastDeletedEvent();

    virtual void getSpecificEventDataString
        (std::stringstream& ss, BranchEvent* event);

    double _updateBetaScale;
    double _updateBetaShiftScale;
    double _updateNodeStateScale;

    double _lastDeletedEventBetaInit;;
    double _lastDeletedEventBetaShift;

    // Here are several variables that track the previous
    // state. At some point, these should have their own class

    // Ultimately, initializations should be handled in TraitModel
    // and Model classes
    // NOT in class TREE!!
    void   initializeTraitParamsForNodes();
    double _lastLH;

    double _readBetaInit;
    double _readBetaShift;
};


inline double TraitModel::getLastLH()
{
    return _lastLH;
}


#endif
