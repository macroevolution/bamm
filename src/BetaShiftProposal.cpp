#include "BetaShiftProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "TraitBranchEvent.h"


BetaShiftProposal::BetaShiftProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(random, settings, model, prior)
{
    _weight = _settings.get<double>("updateRateBetaShift");
    _updateBetaShiftScale = _settings.get<double>("updateBetaShiftScale");
}


double BetaShiftProposal::acceptanceRatio()
{
    if (static_cast<TraitBranchEvent*>(_event)->isTimeVariable()) {
        return EventParameterProposal::acceptanceRatio();
    } else {
        return 0.0;
    }
}


double BetaShiftProposal::getCurrentParameterValue()
{
    return static_cast<TraitBranchEvent*>(_event)->getBetaShift();
}


double BetaShiftProposal::computeNewParameterValue()
{
    return _currentParameterValue + _random.normal(0.0, _updateBetaShiftScale);
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
