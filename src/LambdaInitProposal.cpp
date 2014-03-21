#include "LambdaInitProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "SpExBranchEvent.h"
#include "Tree.h"


LambdaInitProposal::LambdaInitProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        Proposal(rng, settings, model), _prior(prior)
{
    _updateLambdaInitScale = settings.getUpdateLambdaInitScale();
}


void LambdaInitProposal::saveCurrentState()
{
    _event = static_cast<SpExBranchEvent*>(_model.chooseEventAtRandom(true));

    _currentLambdaInit = _event->getLamInit();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
}


void LambdaInitProposal::proposeNewState()
{
    _cterm = std::exp(_updateLambdaInitScale * (_rng.uniformRv() - 0.5));
    _proposedLambdaInit = _cterm * _currentLambdaInit;
    _event->setLamInit(_proposedLambdaInit);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();

    _proposedLogLikelihood = _model.computeLogLikelihood();
}


double LambdaInitProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double LambdaInitProposal::computeLogPriorRatio()
{
    double logPriorRatio = 0.0;

    if (_event == _model.getRootEvent()) {
        logPriorRatio = _prior.lambdaInitRootPrior(_proposedLambdaInit) -
            _prior.lambdaInitRootPrior(_currentLambdaInit);
    } else {
        logPriorRatio = _prior.lambdaInitPrior(_proposedLambdaInit) -
            _prior.lambdaInitPrior(_currentLambdaInit);
    }

    return logPriorRatio;
}


double LambdaInitProposal::computeLogQRatio()
{
    return std::log(_cterm);
}


void LambdaInitProposal::specificAccept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void LambdaInitProposal::specificReject()
{
    _event->setLamInit(_currentLambdaInit);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();
}
