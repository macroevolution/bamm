#include "EventNumberForBranchProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"

#include <algorithm>


EventNumberForBranchProposal::EventNumberForBranchProposal
    (Random& random, Settings& settings, Model& model) :
        _random(random), _model(model)
{
    _weight = settings.get<double>("updateRateEventNumberForBranch");

    _validateEventConfiguration =
        settings.get<bool>("validateEventConfiguration");

    Tree& tree = *(model.getTreePtr());
    _numberOfBranches = tree.getNumberOfNodes() - 1;    // don't count the root
    _totalTreeLength = tree.getTreeLength();
}


void EventNumberForBranchProposal::propose()
{
    _currentEventCount = _model.getNumberOfEvents();
    _currentLogLikelihood = _model.getCurrentLogLikelihood();
    _currentLogPrior = _model.computeLogPrior();

    bool shouldAddEvent = (_currentEventCount == 0) ||
        _random.trueWithProbability(0.5);

    if (shouldAddEvent) {
        _lastEventChanged = _model.addRandomEventToTreeOnRandomBranch();
        _lastProposal = AddEvent;
    } else {
        _lastEventChanged = _model.removeRandomEventFromTree();
        _lastProposal = RemoveEvent;
    }

    _model.setMeanBranchParameters();

    _proposedEventCount = _model.getNumberOfEvents();
    _proposedLogLikelihood = _model.computeLogLikelihood();
    _proposedLogPrior = _model.computeLogPrior();
}


void EventNumberForBranchProposal::accept()
{
    if (_lastProposal == RemoveEvent) {
        if (_lastEventChanged != NULL) {
            delete _lastEventChanged;
            _lastEventChanged = NULL;
        }
    }

    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void EventNumberForBranchProposal::reject()
{
    if (_lastProposal == AddEvent) {
        _model.removeEventFromTree(_lastEventChanged);
        _model.setMeanBranchParameters();
        delete _lastEventChanged;
        _lastEventChanged = NULL;
    } else if (_lastProposal == RemoveEvent) {
        _model.addEventToTree(_lastEventChanged);
    }
}


double EventNumberForBranchProposal::acceptanceRatio()
{
    if (_validateEventConfiguration && _lastProposal == AddEvent &&
            !_model.isEventConfigurationValid(_lastEventChanged)) {
        return 0.0;
    }

    double logLikelihoodRatio = computeLogLikelihoodRatio();
    double logPriorRatio = computeLogPriorRatio();
    double logQRatio = computeLogQRatio();

    double t = _model.getTemperatureMH();
    double logRatio = t * (logLikelihoodRatio + logPriorRatio) + logQRatio;

    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }
}


double EventNumberForBranchProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double EventNumberForBranchProposal::computeLogPriorRatio()
{
    if (_lastProposal == AddEvent) {
        return _proposedLogPrior - _currentLogPrior +
            std::log(_model.getEventRate()) - std::log(_currentEventCount + 1);
    } else {
        return _proposedLogPrior - _currentLogPrior +
            std::log(_currentEventCount) - std::log(_model.getEventRate());
    }
}


double EventNumberForBranchProposal::computeLogQRatio()
{
    double branchLength = _lastEventChanged->getEventNode()->getBrlen();
    double logQRatioBranch =
        std::log((branchLength * _numberOfBranches) / _totalTreeLength);

    if (_lastProposal == AddEvent) {
        // -0.6931... is ln 0.5
        double logQRatio = (_currentEventCount > 0) ? 0.0 : -0.69314718055995;
        return logQRatio - _model.logQRatioJump() + logQRatioBranch;
    } else {
        // 0.6931... is ln 2.0
        double logQRatio = (_currentEventCount != 1) ? 0.0 : 0.69314718055995;
        return logQRatio + _model.logQRatioJump() - logQRatioBranch;
    }
}
