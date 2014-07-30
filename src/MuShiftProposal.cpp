#include "MuShiftProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "SpExBranchEvent.h"


MuShiftProposal::MuShiftProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(random, settings, model, prior)
{
    _weight = _settings.get<double>("updateRateMuShift");
    _updateMuShiftScale = _settings.get<double>("updateMuShiftScale");
}


double MuShiftProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getMuShift();
}


double MuShiftProposal::computeNewParameterValue()
{
    return _currentParameterValue + _random.normal(0.0, _updateMuShiftScale);
}


void MuShiftProposal::setProposedParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setMuShift(_proposedParameterValue);
}


void MuShiftProposal::revertToOldParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setMuShift(_currentParameterValue);
}


void MuShiftProposal::updateParameterOnTree()
{
    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();
}


double MuShiftProposal::computeRootLogPriorRatio()
{
    return _prior.muShiftRootPrior(_proposedParameterValue) -
           _prior.muShiftRootPrior(_currentParameterValue);
}


double MuShiftProposal::computeNonRootLogPriorRatio()
{
    return _prior.muShiftPrior(_proposedParameterValue) -
           _prior.muShiftPrior(_currentParameterValue);
}
