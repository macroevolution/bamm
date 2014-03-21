#include "EventRateProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"


EventRateProposal::EventRateProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        Proposal(rng, settings, model), _prior(prior)
{
    _updateEventRateScale = _settings.getUpdateEventRateScale();
}


void EventRateProposal::saveCurrentState()
{
    _currentEventRate = _model.getEventRate();
}


void EventRateProposal::proposeNewState()
{
    _cterm = std::exp(_updateEventRateScale * (_rng.uniformRv() - 0.5));
    _proposedEventRate = _cterm * _currentEventRate;
    _model.setEventRate(_proposedEventRate);
}


double EventRateProposal::computeLogLikelihoodRatio()
{
    return 0.0;
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


void EventRateProposal::specificAccept()
{
}


void EventRateProposal::specificReject()
{
    _model.setEventRate(_currentEventRate);
}
