#include "Prior.h"
#include "MbRandom.h"
#include "Settings.h"

#define _UPDATE_TOL 0.0001


Prior::Prior(MbRandom* ran, Settings* sp)
{
    ranPtr = ran;
    sttings = sp;
    
    std::cout << "Initialized class Prior object\n\n" << std::endl;
}


Prior::~Prior()
{
}


double Prior::lambdaShiftPrior(double x)
{
    return ranPtr->lnNormalPdf((double)0.0, sttings->getLambdaShiftPrior(), x);
}


double Prior::generateLambdaShiftFromPrior()
{
    if (sttings->getUpdateRateLambdaShift() <= _UPDATE_TOL) {
        return sttings->getLambdaShift0();
    } else {
        return ranPtr->normalRv((double)0.0, sttings->getLambdaShiftPrior());
    }
}


double Prior::lambdaInitPrior(double x)
{
    return ranPtr->lnExponentialPdf(sttings->getLambdaInitPrior(), x);
}


double Prior::generateLambdaInitFromPrior()
{
    if (sttings->getUpdateRateLambda0() <= _UPDATE_TOL) {
        return sttings->getLambdaInit0();
    } else {
        return ranPtr->exponentialRv(sttings->getLambdaInitPrior());
    }
}


double Prior::muInitPrior(double x)
{
    return ranPtr->lnExponentialPdf(sttings->getMuInitPrior(), x);
}


double Prior::generateMuInitFromPrior()
{
    if (sttings->getUpdateRateMu0() <= _UPDATE_TOL){
        return sttings->getMuInit0();
    } else {
        return ranPtr->exponentialRv(sttings->getMuInitPrior());
    }
}


double Prior::muShiftPrior(double x)
{
    return ranPtr->lnNormalPdf((double)0.0, sttings->getMuShiftPrior(), x);
}


double Prior::generateMuShiftFromPrior()
{
    if (sttings->getUpdateRateMuShift() <= _UPDATE_TOL) {
        return sttings->getMuShift0();
    } else {
        return ranPtr->normalRv((double)0.0, sttings->getMuShiftPrior());
    }
}


double Prior::poissonRatePrior(double x)
{
    return ranPtr->lnExponentialPdf(sttings->getPoissonRatePrior(), x);
}


double Prior::generatePoissonRateFromPrior()
{
    return ranPtr->exponentialRv(sttings->getPoissonRatePrior());
}


double Prior::betaInitPrior(double x)
{
    return ranPtr->lnExponentialPdf(sttings->getBetaInitPrior(), x);
}


double Prior::generateBetaInitFromPrior()
{
    if (sttings->getUpdateRateBeta0() <= _UPDATE_TOL){
        return sttings->getBetaInit();
    } else {
        return ranPtr->exponentialRv(sttings->getBetaInitPrior());
    }
}


double Prior::betaShiftPrior(double x)
{
    return ranPtr->lnNormalPdf((double)0.0, sttings->getBetaShiftPrior(), x);
}


double Prior::generateBetaShiftFromPrior()
{
    if (sttings->getUpdateRateBetaShift() <= _UPDATE_TOL){
        return sttings->getBetaShiftInit();
    } else {
        return ranPtr->normalRv((double)0.0, sttings->getBetaShiftPrior());
    }
}
