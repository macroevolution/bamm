#include "MuInitProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "SpExBranchEvent.h"


MuInitProposal::MuInitProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(random, settings, model, prior)
{
    _weight = _settings.get<double>("updateRateMu0");
    _updateMuInitScale = _settings.get<double>("updateMuInitScale");
}


double MuInitProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getMuInit();
}


double MuInitProposal::computeNewParameterValue()
{
    _cterm = std::exp(_updateMuInitScale * (_random.uniform() - 0.5));
    return _cterm * _currentParameterValue;
}


void MuInitProposal::setProposedParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setMuInit(_proposedParameterValue);
}


void MuInitProposal::revertToOldParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setMuInit(_currentParameterValue);
}


void MuInitProposal::updateParameterOnTree()
{
    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();
}


double MuInitProposal::computeRootLogPriorRatio()
{
    return _prior.muInitRootPrior(_proposedParameterValue) -
           _prior.muInitRootPrior(_currentParameterValue);
}


double MuInitProposal::computeNonRootLogPriorRatio()
{
    return _prior.muInitPrior(_proposedParameterValue) -
           _prior.muInitPrior(_currentParameterValue);
}


double MuInitProposal::computeLogQRatio()
{
    return std::log(_cterm);
}
