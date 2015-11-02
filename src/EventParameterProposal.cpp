#include "EventParameterProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "Tree.h"

#include <algorithm>
#include <cmath>


EventParameterProposal::EventParameterProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        _random(random), _settings(settings), _model(model), _prior(prior),
        _tree(_model.getTreePtr())
{
}


void EventParameterProposal::propose()
{
    _event = _model.chooseEventAtRandom(true);
    _currentParameterValue = getCurrentParameterValue();

    _currentLogLikelihood = _model.getCurrentLogLikelihood();

    _proposedParameterValue = computeNewParameterValue();
    setProposedParameterValue();

    updateParameterOnTree();

    _proposedLogLikelihood = _model.computeLogLikelihood();
}


void EventParameterProposal::accept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void EventParameterProposal::reject()
{
    revertToOldParameterValue();
    updateParameterOnTree();
}


double EventParameterProposal::acceptanceRatio()
{
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


double EventParameterProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double EventParameterProposal::computeLogPriorRatio()
{
    if (_event == _model.getRootEvent()) {
        return computeRootLogPriorRatio();
    } else {
        return computeNonRootLogPriorRatio();
    }
}


double EventParameterProposal::computeRootLogPriorRatio()
{
    return 0.0;
}


double EventParameterProposal::computeNonRootLogPriorRatio()
{
    return 0.0;
}


double EventParameterProposal::computeLogQRatio()
{
    return 0.0;
}
