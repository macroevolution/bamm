#include "PreservationRateProposal.h"
#include "Random.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "SpExModel.h"
#include <algorithm>

PreservationRateProposal::PreservationRateProposal
    (Random& random, Settings& settings, Model& model, Prior& prior) :
        _random(random), _settings(settings), _model(model), _prior(prior)
{
    _weight = settings.get<double>("updateRatePreservationRate");
    _updatePreservationRateScale = settings.get<double>("updatePreservationRateScale");
    
}

void PreservationRateProposal::propose()
{
    _currentParameterValue = getCurrentParameterValue();
    
    _currentLogLikelihood = _model.getCurrentLogLikelihood();

    _cterm = std::exp(_updatePreservationRateScale * ( _random.uniform() - 0.5 ));
    _proposedParameterValue = _cterm * _currentParameterValue;
    
    setProposedParameterValue();
 
    _proposedLogLikelihood = _model.computeLogLikelihood();
    
}

double PreservationRateProposal::getCurrentParameterValue()
{
    return static_cast<SpExModel*>(&_model)->getPreservationRate();
}

void PreservationRateProposal::setProposedParameterValue()
{
    static_cast<SpExModel*>(&_model)->setPreservationRate(_proposedParameterValue);
}

void PreservationRateProposal::revertToOldParameterValue()
{
     static_cast<SpExModel*>(&_model)->setPreservationRate(_currentParameterValue);
}


void PreservationRateProposal::accept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void PreservationRateProposal::reject()
{
    revertToOldParameterValue();
}


double PreservationRateProposal::acceptanceRatio()
{
    double logLikelihoodRatio = _proposedLogLikelihood - _currentLogLikelihood;
    double logPriorRatio = computeLogPriorRatio();
    double logQratio = computeLogQRatio();
    
    double t = _model.getTemperatureMH();
    double logRatio = t * (logLikelihoodRatio + logPriorRatio) + logQratio;
    
    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }

}

double PreservationRateProposal::computeLogPriorRatio()
{

    return _prior.preservationRatePrior(_proposedParameterValue) -
                _prior.preservationRatePrior(_currentParameterValue);
    
}


double PreservationRateProposal::computeLogQRatio()
{
    return 0.0;
}




