#include "EventRateProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"

#include <algorithm>

#define USE_ANALYTICAL_POSTERIOR


EventRateProposal::EventRateProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        _random(random), _settings(settings), _model(model), _prior(prior)
{
    _weight = _settings.get<double>("updateRateEventRate");
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
    
    double logQRatio = computeLogQRatio();

//#ifdef USE_ANALYTICAL_POSTERIOR
   
    double logPosteriorRatio = computeLogPosteriorRatio();
    double t = _model.getTemperatureMH();
    double logRatio = t * logPosteriorRatio + logQRatio;

//#else
//    double logPriorRatio = computeLogPriorRatio();   
//    double t = _model.getTemperatureMH();
//    double logRatio = t * logPriorRatio + logQRatio;
//    
//#endif
    
    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }
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


// May 26 , 2015: directly compute log-posterior ratio
//  after correspondence with Bret Larget / Cecile Ane
double EventRateProposal::computeLogPosteriorRatio()
{
    
    double NN = (double)_model.getNumberOfEvents();
    
    double logPosteriorRatio = NN * (std::log(_proposedEventRate) - std::log(_currentEventRate) );
    logPosteriorRatio += (_settings.get<double>("poissonRatePrior") + 1 ) * (_currentEventRate - _proposedEventRate);
    
    return logPosteriorRatio;

}

#undef USE_ANALYTICAL_POSTERIOR


