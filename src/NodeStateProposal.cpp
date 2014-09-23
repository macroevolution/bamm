#include "NodeStateProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "TraitModel.h"
#include "Tree.h"
#include "Node.h"
#include "Stat.h"

#include <cstdlib>
#include <algorithm>


NodeStateProposal::NodeStateProposal
    (Random& random, Settings& settings, Model& model) :
        _random(random), _settings(settings),
        _model(static_cast<TraitModel&>(model)), _tree(model.getTreePtr())
{
    _weight = _settings.get<double>("updateRateNodeState");

    // Node state scale is relative to the standard deviation
    // of the trait values (located in the tree terminal nodes)
    double sd_traits = Stat::standard_deviation(_tree->traitValues());
    _updateNodeStateScale =
        _settings.get<double>("updateNodeStateScale") * sd_traits;

    _priorMin = _settings.get<double>("traitPriorMin");
    _priorMax = _settings.get<double>("traitPriorMax");

    // Min and max trait priors must be updated later,
    // after tree is initialized during model construction
    _minMaxTraitPriorUpdated = false;
}


void NodeStateProposal::updateMinMaxTraitPriorSettings()
{
    int nnodes = _tree->getNumberOfNodes();
    std::vector<double> tvec;

    const std::vector<Node*>& postOrderNodes = _tree->postOrderNodes();
    for (int i = 0; i < nnodes; i++) {
        Node* xnode = postOrderNodes[i];
        if (xnode->getTraitValue() != 0) {
            tvec.push_back(xnode->getTraitValue());
        }
    }

    std::sort(tvec.begin(), tvec.end());

    // Default here will be to use observed range +/- 20%
    double rg = tvec[(tvec.size() - 1)] - tvec[0];
    _priorMin = tvec[0] - (0.2 * rg);
    _priorMax = tvec[(tvec.size() - 1)] + (0.2 * rg);
}


void NodeStateProposal::propose()
{
    if (!_minMaxTraitPriorUpdated) {
        updateMinMaxTraitPriorSettings();
        _minMaxTraitPriorUpdated = true;
    }

    _node = _tree->chooseInternalNodeAtRandom();

    double currentTriadLogLikelihood =
        _model.computeTriadLikelihoodTraits(_node);
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
    _currentNodeState = _node->getTraitValue();

    _proposedNodeState = _currentNodeState + _random.uniform
        (-_updateNodeStateScale, _updateNodeStateScale);
    _node->setTraitValue(_proposedNodeState);

    double proposedTriadLogLikelihood =
        _model.computeTriadLikelihoodTraits(_node);
    _proposedLogLikelihood = _currentLogLikelihood -
        currentTriadLogLikelihood + proposedTriadLogLikelihood;
}


void NodeStateProposal::accept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void NodeStateProposal::reject()
{
    _node->setTraitValue(_currentNodeState);
}


double NodeStateProposal::acceptanceRatio()
{
    if (_proposedNodeState < _priorMin || _proposedNodeState > _priorMax) {
        return 0.0;
    }

    double logLikelihoodRatio = computeLogLikelihoodRatio();

    double t = _model.getTemperatureMH();
    double logRatio = t * logLikelihoodRatio;

    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }
}


double NodeStateProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}
