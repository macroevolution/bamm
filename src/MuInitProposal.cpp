#include "MuInitProposal.h"
#include "MbRandom.h"
#include "Settings.h"
#include "Model.h"
#include "Prior.h"
#include "SpExBranchEvent.h"


MuInitProposal::MuInitProposal
    (MbRandom& rng, Settings& settings, Model& model, Prior& prior) :
        EventParameterProposal(rng, settings, model, prior)
{
    _updateMuInitScale = _settings.getUpdateMuInitScale();
}


double MuInitProposal::getCurrentParameterValue()
{
    return static_cast<SpExBranchEvent*>(_event)->getMuInit();
}


double MuInitProposal::computeNewParameterValue()
{
    _cterm = std::exp(_updateMuInitScale * (_rng.uniformRv() - 0.5));
    return _cterm * _currentParameterValue;
}


void MuInitProposal::setProposedParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setMuInit(_proposedParameterValue);
}


void MuInitProposal::revertToOldParameterValue()
{
    static_cast<SpExBranchEvent*>(_event)->setMuInit(_currentParameterValue);
}


double MuInitProposal::computeRootLogPriorRatio()
{
    return _prior.muInitRootPrior(_proposedParameterValue) -
           _prior.muInitRootPrior(_currentParameterValue);
}


double MuInitProposal::computeNonRootLogPriorRatio()
{
    return _prior.muInitPrior(_proposedParameterValue) -
           _prior.muInitPrior(_currentParameterValue);
}


double MuInitProposal::computeLogQRatio()
{
    return std::log(_cterm);
}
