#include "MoveEventProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Node.h"
#include "BranchHistory.h"
#include "Tree.h"

#include <algorithm>


MoveEventProposal::MoveEventProposal
    (Random& random, Settings& settings, Model& model) :
        _random(random), _settings(settings), _model(model)
{
    _weight = _settings.get<double>("updateRateEventPosition");

    _localToGlobalMoveRatio = _settings.get<double>("localGlobalMoveRatio");
    _scale = _settings.get<double>("updateEventLocationScale") *
        _model.getTreePtr()->maxRootToTipLength();

    _validateEventConfiguration =
        _settings.get<bool>("validateEventConfiguration");
}


void MoveEventProposal::propose()
{
    _currentEventCount = _model.getNumberOfEvents();
    if (_currentEventCount == 0) {
        return;
    }

    _event = _model.chooseEventAtRandom();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();

    // This is the event preceding the chosen event;
    // histories should be set forward from here
    BranchEvent* previousEvent = _event->getEventNode()->getBranchHistory()->
        getLastEvent(_event);

    _event->getEventNode()->getBranchHistory()->
        popEventOffBranchHistory(_event);

    double localMoveProb = _localToGlobalMoveRatio /
        (1 + _localToGlobalMoveRatio);

    // Choose to move locally or globally
    if (_random.trueWithProbability(localMoveProb)) {
        double step = _random.uniform(0, _scale) - 0.5 * _scale;
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


void MoveEventProposal::accept()
{
    if (_currentEventCount == 0) {
        return;
    }

    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void MoveEventProposal::reject()
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


double MoveEventProposal::acceptanceRatio()
{
    if (_currentEventCount == 0) {
        return 0.0;
    }

    if (_validateEventConfiguration &&
            !_model.isEventConfigurationValid(_event)) {
        return 0.0;
    }

    double logLikelihoodRatio = computeLogLikelihoodRatio();

    double t = _model.getTemperatureMH();
    double logRatio = t * logLikelihoodRatio;

    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }
}


double MoveEventProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}
