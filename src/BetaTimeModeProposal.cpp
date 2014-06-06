#include "BetaTimeModeProposal.h"

#include "Prior.h"
#include "TraitBranchEvent.h"

class Random;
class Settings;
class Model;


BetaTimeModeProposal::BetaTimeModeProposal
    (Random& random, Settings& settings, Model& model)
    : TimeModeProposal(random, settings, model)
{
}


double BetaTimeModeProposal::initialParameter(BranchEvent* event)
{
    return static_cast<TraitBranchEvent*>(event)->getBetaInit();
}


double BetaTimeModeProposal::rateParameter(BranchEvent* event)
{
    return static_cast<TraitBranchEvent*>(event)->getBetaShift();
}


bool BetaTimeModeProposal::isTimeVariable(BranchEvent* event)
{
    return static_cast<TraitBranchEvent*>(event)->isTimeVariable();
}


void BetaTimeModeProposal::setEventParameters(BranchEvent* event,
    double initParam, double rateParam, bool isTimeVariable)
{
    TraitBranchEvent* traitEvent = static_cast<TraitBranchEvent*>(event);

    traitEvent->setBetaInit(initParam);
    traitEvent->setBetaShift(rateParam);
    traitEvent->setTimeVariable(isTimeVariable);
}


void BetaTimeModeProposal::setModelParameters()
{
    _tree->setMeanBranchTraitRates();
}


double BetaTimeModeProposal::rateParameterFromPrior()
{
    return _prior.generateBetaShiftFromPrior();
}
