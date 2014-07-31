#include "LambdaShiftProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "SpExBranchEvent.h"


LambdaShiftProposal::LambdaShiftProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(random, settings, model, prior)
{
    _weight = _settings.get<double>("updateRateLambdaShift");
    _updateLambdaShiftScale = _settings.get<double>("updateLambdaShiftScale");
}


double LambdaShiftProposal::acceptanceRatio()
{
    if (static_cast<SpExBranchEvent*>(_event)->isTimeVariable()) {
        return EventParameterProposal::acceptanceRatio();
    } else {
        return 0.0;
    }
}


double LambdaShiftProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getLamShift();
}


double LambdaShiftProposal::computeNewParameterValue()
{
    return _currentParameterValue +
        _random.normal(0.0, _updateLambdaShiftScale);
}


void LambdaShiftProposal::setProposedParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setLamShift(_proposedParameterValue);
}


void LambdaShiftProposal::revertToOldParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setLamShift(_currentParameterValue);
}


void LambdaShiftProposal::updateParameterOnTree()
{
    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();
}


double LambdaShiftProposal::computeRootLogPriorRatio()
{
    return _prior.lambdaShiftRootPrior(_proposedParameterValue) -
           _prior.lambdaShiftRootPrior(_currentParameterValue);
}


double LambdaShiftProposal::computeNonRootLogPriorRatio()
{
    return _prior.lambdaShiftPrior(_proposedParameterValue) -
           _prior.lambdaShiftPrior(_currentParameterValue);
}
