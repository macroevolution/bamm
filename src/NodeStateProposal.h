#ifndef NODE_STATE_PROPOSAL_H
#define NODE_STATE_PROPOSAL_H


#include "Proposal.h"

class Random;
class Settings;
class Model;
class TraitModel;
class Tree;
class Node;


class NodeStateProposal : public Proposal
{
public:

    NodeStateProposal(Random& random, Settings& settings, Model& model);

    virtual void propose();
    virtual void accept();
    virtual void reject();

    virtual double acceptanceRatio();

private:

    void updateMinMaxTraitPriorSettings();
    double computeLogLikelihoodRatio();

    Random& _random;
    Settings& _settings;
    TraitModel& _model;

    Tree* _tree;
    Node* _node;

    double _updateNodeStateScale;
    double _priorMin;
    double _priorMax;

    bool _minMaxTraitPriorUpdated;

    double _currentNodeState;
    double _proposedNodeState;

    double _currentLogLikelihood;
    double _proposedLogLikelihood;
};


#endif
