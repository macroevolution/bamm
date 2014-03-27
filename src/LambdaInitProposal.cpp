#include "LambdaInitProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "SpExBranchEvent.h"


LambdaInitProposal::LambdaInitProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(rng, settings, model, prior)
{
    _updateLambdaInitScale = _settings.getUpdateLambdaInitScale();
}


double LambdaInitProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getLamInit();
}


double LambdaInitProposal::computeNewParameterValue()
{
    _cterm = std::exp(_updateLambdaInitScale * (_rng.uniformRv() - 0.5));
    return _cterm * _currentParameterValue;
}


void LambdaInitProposal::setProposedParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setLamInit(_proposedParameterValue);
}


void LambdaInitProposal::revertToOldParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setLamInit(_currentParameterValue);
}


double LambdaInitProposal::computeRootLogPriorRatio()
{
    return _prior.lambdaInitRootPrior(_proposedParameterValue) -
           _prior.lambdaInitRootPrior(_currentParameterValue);
}


double LambdaInitProposal::computeNonRootLogPriorRatio()
{
    return _prior.lambdaInitPrior(_proposedParameterValue) -
           _prior.lambdaInitPrior(_currentParameterValue);
}


double LambdaInitProposal::computeLogQRatio()
{
    return std::log(_cterm);
}
