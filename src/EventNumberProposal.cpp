#include "EventNumberProposal.h"
#include "MbRandom.h"
#include "Model.h"
#include "Node.h"
#include "BranchHistory.h"


EventNumberProposal::EventNumberProposal
    (MbRandom& rng, Settings& settings, Model& model) :
        Proposal(rng, settings, model)
{
}


void EventNumberProposal::saveCurrentState()
{
    _currentEventCount = _model.getNumberOfEvents();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
    _currentLogPrior = _model.computeLogPrior();
}


void EventNumberProposal::proposeNewState()
{
    bool shouldAddEvent = _rng.uniformRv() < 0.5;
    if (shouldAddEvent || _model.getNumberOfEvents() == 0) {
        _lastProposal = AddEvent;
        proposeAddEvent();
    } else {
        _lastProposal = RemoveEvent;
        proposeRemoveEvent();
    }
}


void EventNumberProposal::proposeAddEvent()
{
    _lastEventChanged = _model.addRandomEventToTree();
    _model.setMeanBranchParameters();

    _proposedEventCount = _currentEventCount + 1;
    _proposedLogLikelihood = _model.computeLogLikelihood();
    _proposedLogPrior = _model.computeLogPrior();
}


void EventNumberProposal::proposeRemoveEvent()
{
    // Flag proposal for rejection if there are no events to remove
    if (_currentEventCount == 0) {
        _model.setProposalFail(true);
    } else {
        _model.setProposalFail(false);
    }

    _lastEventChanged = _model.removeRandomEventFromTree();
    _model.setMeanBranchParameters();

    _proposedEventCount = _currentEventCount - 1;
    _proposedLogLikelihood = _model.computeLogLikelihood();
    _proposedLogPrior = _model.computeLogPrior();
}


double EventNumberProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double EventNumberProposal::computeLogPriorRatio()
{
    if (_lastProposal == AddEvent) {
        return _proposedLogPrior - _currentLogPrior +
            std::log(_model.getEventRate()) - std::log(_proposedEventCount);
    } else {
        return _proposedLogPrior - _currentLogPrior +
            std::log(_proposedEventCount) - std::log(_model.getEventRate());
    }
}


double EventNumberProposal::computeLogQRatio()
{
    if (_lastProposal == AddEvent) {
        // -0.6931... is ln 0.5
        double logQRatio = (_currentEventCount > 0) ? 0.0 : -0.69314718055995;
        return logQRatio - _model.logQRatioJump();
    } else {
        // 0.6931... is ln 2.0
        double logQRatio = (_currentEventCount != 1) ? 0.0 : 0.69314718055995;
        return logQRatio + _model.logQRatioJump();
    }
}


void EventNumberProposal::specificAccept()
{
    if (_lastProposal == RemoveEvent) {
        if (_lastEventChanged != NULL) {
            delete _lastEventChanged;
            _lastEventChanged = NULL;
        }
    }

    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void EventNumberProposal::specificReject()
{
    if (_lastProposal == AddEvent) {
        rejectAddEvent();
    } else if (_lastProposal == RemoveEvent) {
        rejectRemoveEvent();
    }
}


void EventNumberProposal::rejectAddEvent()
{
    _model.removeEventFromTree(_lastEventChanged);

    delete _lastEventChanged;
    _lastEventChanged = NULL;

    _model.setMeanBranchParameters();
}


void EventNumberProposal::rejectRemoveEvent()
{
    _model.addEventToTree(_lastEventChanged);

    // Add event to branch history
    _lastEventChanged->getEventNode()->getBranchHistory()->
        addEventToBranchHistory(_lastEventChanged);

    _model.events().insert(_lastEventChanged);

    // Event is inserted into branch history,
    // but branch histories must be updated
    _model.forwardSetBranchHistories(_lastEventChanged);

    _model.setMeanBranchParameters();
}
