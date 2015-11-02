#include "BetaInitProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"
#include "TraitBranchEvent.h"


BetaInitProposal::BetaInitProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(random, settings, model, prior)
{
    _weight = _settings.get<double>("updateRateBeta0");
    _updateBetaInitScale = _settings.get<double>("updateBetaInitScale");
}


double BetaInitProposal::getCurrentParameterValue()
{
    return static_cast<TraitBranchEvent*>(_event)->getBetaInit();
}


double BetaInitProposal::computeNewParameterValue()
{
    _cterm = std::exp(_updateBetaInitScale * (_random.uniform() - 0.5));
    return _cterm * _currentParameterValue;
}


void BetaInitProposal::setProposedParameterValue()
{
    static_cast<TraitBranchEvent*>(_event)->
        setBetaInit(_proposedParameterValue);
}


void BetaInitProposal::revertToOldParameterValue()
{
    static_cast<TraitBranchEvent*>(_event)->setBetaInit(_currentParameterValue);
}


void BetaInitProposal::updateParameterOnTree()
{
    _tree->setMeanBranchTraitRates();
}


double BetaInitProposal::computeRootLogPriorRatio()
{
    return _prior.betaInitRootPrior(_proposedParameterValue) -
           _prior.betaInitRootPrior(_currentParameterValue);
}


double BetaInitProposal::computeNonRootLogPriorRatio()
{
    return _prior.betaInitPrior(_proposedParameterValue) -
           _prior.betaInitPrior(_currentParameterValue);
}


double BetaInitProposal::computeLogQRatio()
{
    return std::log(_cterm);
}
