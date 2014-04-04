#include "NodeStateProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "TraitModel.h"
#include "Tree.h"
#include "Node.h"
#include "Stat.h"

#include <cstdlib>
#include <algorithm>


NodeStateProposal::NodeStateProposal
    (MbRandom& rng, Settings& settings, Model& model) :
        _rng(rng), _settings(settings), _model(static_cast<TraitModel&>(model)),
        _tree(model.getTreePtr())
{
    // Node state scale is relative to the standard deviation
    // of the trait values (located in the tree terminal nodes)
    double sd_traits = Stat::standard_deviation(_tree->traitValues());
    _updateNodeStateScale = _settings.getUpdateNodeStateScale() * sd_traits;

    updateMinMaxTraitPriorSettings();

    _priorMin = _settings.getTraitPriorMin();
    _priorMax = _settings.getTraitPriorMax();
}


void NodeStateProposal::updateMinMaxTraitPriorSettings()
{
    int nnodes = _tree->getNumberOfNodes();
    std::vector<double> tvec;
    for (int i = 0; i < nnodes; i++) {
        Node* xnode = _tree->getNodeFromDownpassSeq(i);
        if (xnode->getTraitValue() != 0) {
            tvec.push_back(xnode->getTraitValue());
        }
    }

    std::sort(tvec.begin(), tvec.end());

    // Default here will be to use observed range +/- 20%
    double rg = tvec[(tvec.size() - 1)] - tvec[0];
    double minprior = tvec[0] - (0.2 * rg);
    double maxprior = tvec[(tvec.size() - 1)] + (0.2 * rg);

    log() << "\nMin and max phenotype limits set using observed data:\n"
          << "\t\tMin: " << minprior << "\tMax: " << maxprior << "\n";

    _settings.setTraitPriorMin(minprior);
    _settings.setTraitPriorMax(maxprior);
}


void NodeStateProposal::propose()
{
    _node = _tree->chooseInternalNodeAtRandom();

    double currentTriadLogLikelihood =
        _model.computeTriadLikelihoodTraits(_node);
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
    _currentNodeState = _node->getTraitValue();

    _proposedNodeState = _currentNodeState + _rng.uniformRv
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

    return std::min(1.0, std::exp(logRatio));
}


double NodeStateProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}
