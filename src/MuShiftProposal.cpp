#include "MuShiftProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "SpExBranchEvent.h"


MuShiftProposal::MuShiftProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(rng, settings, model, prior)
{
    _updateMuShiftScale = _settings.getUpdateMuShiftScale();
}


double MuShiftProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getMuShift();
}


double MuShiftProposal::computeNewParameterValue()
{
    return _currentParameterValue + _rng.normalRv(0.0, _updateMuShiftScale);
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
