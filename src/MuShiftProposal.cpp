#include "MuShiftProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "SpExBranchEvent.h"
#include "Tree.h"


MuShiftProposal::MuShiftProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        Proposal(rng, settings, model), _prior(prior)
{
    _updateMuShiftScale = settings.getUpdateMuShiftScale();
}


void MuShiftProposal::saveCurrentState()
{
    _event = static_cast<SpExBranchEvent*>(_model.chooseEventAtRandom(true));

    _currentMuShift = _event->getMuShift();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
}


void MuShiftProposal::proposeNewState()
{
    _proposedMuShift = _currentMuShift +
        _rng.normalRv(0.0, _updateMuShiftScale);
    _event->setMuShift(_proposedMuShift);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();

    _proposedLogLikelihood = _model.computeLogLikelihood();
}


double MuShiftProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double MuShiftProposal::computeLogPriorRatio()
{
    double logPriorRatio = 0.0;

    if (_event == _model.getRootEvent()) {
        logPriorRatio = _prior.muShiftRootPrior(_proposedMuShift) -
            _prior.muShiftRootPrior(_currentMuShift);
    } else {
        logPriorRatio = _prior.muShiftPrior(_proposedMuShift) -
            _prior.muShiftPrior(_currentMuShift);
    }

    return logPriorRatio;
}


double MuShiftProposal::computeLogQRatio()
{
    return 0.0;
}


void MuShiftProposal::specificAccept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void MuShiftProposal::specificReject()
{
    _event->setMuShift(_currentMuShift);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();
}
