#include "LambdaShiftProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "SpExBranchEvent.h"


LambdaShiftProposal::LambdaShiftProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(rng, settings, model, prior)
{
    _updateLambdaShiftScale = _settings.getUpdateLambdaShiftScale();
}


double LambdaShiftProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getLamShift();
}


double LambdaShiftProposal::computeNewParameterValue()
{
    return _currentParameterValue + _rng.normalRv(0.0, _updateLambdaShiftScale);
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
