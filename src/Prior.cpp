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
    if (sttings->getUpdateRateMu0() <= _UPDATE_TOL) {
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
    if (sttings->getUpdateRateBeta0() <= _UPDATE_TOL) {
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
    if (sttings->getUpdateRateBetaShift() <= _UPDATE_TOL) {
        return sttings->getBetaShiftInit();
    } else {
        return ranPtr->normalRv((double)0.0, sttings->getBetaShiftPrior());
    }
}


double Prior::lambdaInitRootPrior(double x)
{
    if (fabs(sttings->getLambdaInitRootPrior() + 1) < _UPDATE_TOL) {
        return lambdaInitPrior(x);
    } else {
        return ranPtr->lnExponentialPdf(sttings->getLambdaInitRootPrior(), x);
    }

}


double Prior::lambdaShiftRootPrior(double x)
{
    if (fabs(sttings->getLambdaShiftRootPrior() + 1) < _UPDATE_TOL) {
        return lambdaShiftPrior(x);
    } else {
        return ranPtr->lnExponentialPdf(sttings->getLambdaShiftRootPrior(), x);
    }
}


double Prior::muInitRootPrior(double x)
{
    if (fabs(sttings->getMuInitRootPrior() + 1) < _UPDATE_TOL) {
        return muInitPrior(x);
    } else {
        return ranPtr->lnExponentialPdf(sttings->getMuInitRootPrior(), x);
    }
}


double Prior::muShiftRootPrior(double x)
{
    if (fabs(sttings->getMuShiftRootPrior() + 1) < _UPDATE_TOL) {
        return muShiftPrior(x);
    } else {
        return ranPtr->lnExponentialPdf(sttings->getMuShiftRootPrior(), x);
    }
}


double Prior::betaInitRootPrior(double x)
{
    if (fabs(sttings->getBetaInitRootPrior() + 1) < _UPDATE_TOL) {
        return betaInitPrior(x);
    } else {
        return ranPtr->lnExponentialPdf(sttings->getBetaInitRootPrior(), x);
    }
}


double Prior::betaShiftRootPrior(double x)
{
    if (fabs(sttings->getBetaShiftRootPrior() + 1) < _UPDATE_TOL) {
        return betaShiftPrior(x);
    } else {
        return ranPtr->lnExponentialPdf(sttings->getBetaShiftRootPrior(), x);
    }
}