#include "EventRateProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"


EventRateProposal::EventRateProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        _random(random), _settings(settings), _model(model), _prior(prior)
{
    _updateEventRateScale = _settings.get<double>("updateEventRateScale");
}


void EventRateProposal::propose()
{
    _currentEventRate = _model.getEventRate();

    _cterm = std::exp(_updateEventRateScale * (_random.uniform() - 0.5));
    _proposedEventRate = _cterm * _currentEventRate;

    _model.setEventRate(_proposedEventRate);
}


void EventRateProposal::accept()
{
}


void EventRateProposal::reject()
{
    _model.setEventRate(_currentEventRate);
}


double EventRateProposal::acceptanceRatio()
{
    double logPriorRatio = computeLogPriorRatio();
    double logQRatio = computeLogQRatio();

    double t = _model.getTemperatureMH();
    double logRatio = t * logPriorRatio + logQRatio;

    return std::min(1.0, std::exp(logRatio));
}


double EventRateProposal::computeLogPriorRatio()
{
    return _prior.poissonRatePrior(_proposedEventRate) -
           _prior.poissonRatePrior(_currentEventRate);
}


double EventRateProposal::computeLogQRatio()
{
    return std::log(_cterm);
}
