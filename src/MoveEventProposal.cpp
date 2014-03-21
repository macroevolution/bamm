#include "MoveEventProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Node.h"
#include "BranchHistory.h"
#include "Tree.h"


MoveEventProposal::MoveEventProposal
    (MbRandom& rng, Settings& settings, Model& model) :
        Proposal(rng, settings, model)
{
    _localToGlobalMoveRatio = _settings.getLocalGlobalMoveRatio();
    _scale = _settings.getUpdateEventLocationScale() *
        _model.getTreePtr()->maxRootToTipLength();
}


void MoveEventProposal::saveCurrentState()
{
    _currentEventCount = _model.getNumberOfEvents();
    if (_currentEventCount == 0) {
        return;
    }

    _event = _model.chooseEventAtRandom();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
}


void MoveEventProposal::proposeNewState()
{
    if (_currentEventCount == 0) {
        _model.setProposalFail(true);
        return;
    } else {
        _model.setProposalFail(false);
    }

    // This is the event preceding the chosen event;
    // histories should be set forward from here
    BranchEvent* previousEvent = _event->getEventNode()->getBranchHistory()->
        getLastEvent(_event);

    _event->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_event);

    double localMoveProb = _localToGlobalMoveRatio /
        (1 + _localToGlobalMoveRatio);

    // Choose to move locally or globally
    if (_rng.uniformRv() < localMoveProb) {
        double step = _rng.uniformRv(0, _scale) - 0.5 * _scale;
        _event->moveEventLocal(step);
    } else {
        _event->moveEventGlobal();
    }

    _event->getEventNode()->getBranchHistory()->addEventToBranchHistory(_event);

    _model.forwardSetBranchHistories(previousEvent);
    _model.forwardSetBranchHistories(_event);
    _model.setMeanBranchParameters();

    _proposedLogLikelihood = _model.computeLogLikelihood();
}


double MoveEventProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double MoveEventProposal::computeLogPriorRatio()
{
    return 0.0;
}


double MoveEventProposal::computeLogQRatio()
{
    return 0.0;
}


void MoveEventProposal::specificAccept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void MoveEventProposal::specificReject()
{
    if (_currentEventCount == 0) {
        return;
    }

    // Get last event from position of event to be removed
    BranchEvent* lastEvent = _event->getEventNode()->getBranchHistory()->
        getLastEvent(_event);

    // Pop event off its new location
    _event->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_event);

    // Reset nodeptr, reset mapTime
    _event->revertOldMapPosition();

    // Now reset forward from _lastEventChanged (new position)
    // and from newLastEvent, which holds 'last' event before old position
    _event->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(_event);

    _model.forwardSetBranchHistories(lastEvent);
    _model.forwardSetBranchHistories(_event);
    _model.setMeanBranchParameters();
}
