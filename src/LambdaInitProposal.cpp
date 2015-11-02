#include "LambdaInitProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "SpExBranchEvent.h"


LambdaInitProposal::LambdaInitProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(random, settings, model, prior)
{
    _weight = _settings.get<double>("updateRateLambda0");
    _updateLambdaInitScale = _settings.get<double>("updateLambdaInitScale");
}


double LambdaInitProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getLamInit();
}


double LambdaInitProposal::computeNewParameterValue()
{
    _cterm = std::exp(_updateLambdaInitScale * (_random.uniform() - 0.5));
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


void LambdaInitProposal::updateParameterOnTree()
{
    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();
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
