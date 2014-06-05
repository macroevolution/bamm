#include "LambdaTimeModeProposal.h"

#include "Prior.h"
#include "SpExBranchEvent.h"

class Random;
class Settings;
class Model;


LambdaTimeModeProposal::LambdaTimeModeProposal
    (Random& random, Settings& settings, Model& model)
    : TimeModeProposal(random, settings, model)
{
}


double LambdaTimeModeProposal::initialParameter(BranchEvent* event)
{
    return static_cast<SpExBranchEvent*>(event)->getLamInit();
}


double LambdaTimeModeProposal::rateParameter(BranchEvent* event)
{
    return static_cast<SpExBranchEvent*>(event)->getLamShift();
}


bool LambdaTimeModeProposal::isTimeVariable(BranchEvent* event)
{
    return static_cast<SpExBranchEvent*>(event)->isTimeVariable();
}


void LambdaTimeModeProposal::setEventParameters(BranchEvent* event,
    double initParam, double rateParam, bool isTimeVariable)
{
    SpExBranchEvent* spExEvent = static_cast<SpExBranchEvent*>(event);

    spExEvent->setLamInit(initParam);
    spExEvent->setLamShift(rateParam);
    spExEvent->setTimeVariable(isTimeVariable);
}


double LambdaTimeModeProposal::rateParameterFromPrior()
{
    return _prior.generateLambdaShiftFromPrior();
}
