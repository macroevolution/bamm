#include "Proposal.h"
#include "MbRandom.h"
#include "Model.h"


Proposal::Proposal(MbRandom& rng, Settings& settings, Model& model) :
    _rng(rng), _settings(settings), _model(model)
{
}


Proposal::~Proposal()
{
}


void Proposal::propose()
{
    saveCurrentState();
    proposeNewState();
    updateLogRatios();
}


void Proposal::updateLogRatios()
{
    _model.setLogLikelihoodRatio(computeLogLikelihoodRatio());
    _model.setLogPriorRatio(computeLogPriorRatio());
    _model.setLogQRatio(computeLogQRatio());
}


// TODO: Somewhere I need to check whether event configuration is valid
void Proposal::accept()
{
    specificAccept();
}


void Proposal::reject()
{
    specificReject();
}
