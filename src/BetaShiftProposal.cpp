#include "BetaShiftProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "TraitBranchEvent.h"


BetaShiftProposal::BetaShiftProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(rng, settings, model, prior)
{
    _updateBetaShiftScale = _settings.getUpdateBetaShiftScale();
}


double BetaShiftProposal::getCurrentParameterValue()
{
    return static_cast<TraitBranchEvent*>(_event)->getBetaShift();
}


double BetaShiftProposal::computeNewParameterValue()
{
    return _currentParameterValue + _rng.normalRv(0.0, _updateBetaShiftScale);
}


void BetaShiftProposal::setProposedParameterValue()
{
    static_cast<TraitBranchEvent*>(_event)->
        setBetaShift(_proposedParameterValue);
}


void BetaShiftProposal::revertToOldParameterValue()
{
    static_cast<TraitBranchEvent*>(_event)->
        setBetaShift(_currentParameterValue);
}


void BetaShiftProposal::updateParameterOnTree()
{
    _tree->setMeanBranchTraitRates();
}


double BetaShiftProposal::computeRootLogPriorRatio()
{
    return _prior.betaShiftRootPrior(_proposedParameterValue) -
           _prior.betaShiftRootPrior(_currentParameterValue);
}


double BetaShiftProposal::computeNonRootLogPriorRatio()
{
    return _prior.betaShiftPrior(_proposedParameterValue) -
           _prior.betaShiftPrior(_currentParameterValue);
}
