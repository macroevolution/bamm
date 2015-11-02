#ifndef LAMBDA_TIME_MODE_PROPOSAL_H
#define LAMBDA_TIME_MODE_PROPOSAL_H


#include "TimeModeProposal.h"

class Random;
class Settings;
class Model;


class LambdaTimeModeProposal : public TimeModeProposal
{
public:

    LambdaTimeModeProposal(Random& random, Settings& settings, Model& model);

protected:

    double initialParameter(BranchEvent* event);
    double rateParameter(BranchEvent* event);
    bool isTimeVariable(BranchEvent* event);
    
    void setEventParameters(BranchEvent* event,
        double initParam, double rateParam, bool isTimeVariable);

    void setModelParameters();

    double rateParameterFromPrior();

};


#endif
