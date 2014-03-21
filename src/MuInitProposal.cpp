#include "MuInitProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "SpExBranchEvent.h"
#include "Tree.h"


MuInitProposal::MuInitProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        Proposal(rng, settings, model), _prior(prior)
{
    _updateMuInitScale = settings.getUpdateMuInitScale();
}


void MuInitProposal::saveCurrentState()
{
    _event = static_cast<SpExBranchEvent*>(_model.chooseEventAtRandom(true));

    _currentMuInit = _event->getMuInit();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
}


void MuInitProposal::proposeNewState()
{
    _cterm = std::exp(_updateMuInitScale * (_rng.uniformRv() - 0.5));
    _proposedMuInit = _cterm * _currentMuInit;
    _event->setMuInit(_proposedMuInit);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();

    _proposedLogLikelihood = _model.computeLogLikelihood();
}


double MuInitProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double MuInitProposal::computeLogPriorRatio()
{
    double logPriorRatio = 0.0;

    if (_event == _model.getRootEvent()) {
        logPriorRatio = _prior.muInitRootPrior(_proposedMuInit) -
            _prior.muInitRootPrior(_currentMuInit);
    } else {
        logPriorRatio = _prior.muInitPrior(_proposedMuInit) -
            _prior.muInitPrior(_currentMuInit);
    }

    return logPriorRatio;
}


double MuInitProposal::computeLogQRatio()
{
    return std::log(_cterm);
}


void MuInitProposal::specificAccept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void MuInitProposal::specificReject()
{
    _event->setMuInit(_currentMuInit);

    _model.getTreePtr()->setNodeSpeciationParameters();
    _model.getTreePtr()->setNodeExtinctionParameters();
}
