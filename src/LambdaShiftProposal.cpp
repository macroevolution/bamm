#include "LambdaShiftProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "SpExBranchEvent.h"
#include "Tree.h"


LambdaShiftProposal::LambdaShiftProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        Proposal(rng, settings, model), _prior(prior)
{
    _updateLambdaShiftScale = settings.getUpdateLambdaShiftScale();
}


void LambdaShiftProposal::saveCurrentState()
{
    _event = static_cast<SpExBranchEvent*>(_model.chooseEventAtRandom(true));

    _currentLambdaShift = _event->getLamShift();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
}


void LambdaShiftProposal::proposeNewState()
{
    _proposedLambdaShift = _currentLambdaShift +
        _rng.normalRv(0.0, _updateLambdaShiftScale);
    _event->setLamShift(_proposedLambdaShift);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();

    _proposedLogLikelihood = _model.computeLogLikelihood();
}


double LambdaShiftProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double LambdaShiftProposal::computeLogPriorRatio()
{
    double logPriorRatio = 0.0;

    if (_event == _model.getRootEvent()) {
        logPriorRatio = _prior.lambdaShiftRootPrior(_proposedLambdaShift) -
            _prior.lambdaShiftRootPrior(_currentLambdaShift);
    } else {
        logPriorRatio = _prior.lambdaShiftPrior(_proposedLambdaShift) -
            _prior.lambdaShiftPrior(_currentLambdaShift);
    }

    return logPriorRatio;
}


double LambdaShiftProposal::computeLogQRatio()
{
    return 0.0;
}


void LambdaShiftProposal::specificAccept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void LambdaShiftProposal::specificReject()
{
    _event->setLamShift(_currentLambdaShift);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();
}
