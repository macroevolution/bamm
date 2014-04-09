#include "Prior.h"
#include "MbRandom.h"
#include "Settings.h"

#define _UPDATE_TOL 0.0001


Prior::Prior(MbRandom* rng, Settings* settings) : _rng(rng)
{
    std::string modelType = settings->get("modeltype");

    if (modelType == "speciationextinction") {
        _lambdaInit0 = settings->get<double>("lambdaInit0");
        _muInit0 = settings->get<double>("muInit0");
        _lambdaShift0 = settings->get<double>("lambdaShift0");
        _muShift0 = settings->get<double>("muShift0");
        _lambdaInitPrior = settings->get<double>("lambdaInitPrior");
        _muInitPrior = settings->get<double>("muInitPrior");
        _lambdaShiftPrior = settings->get<double>("lambdaShiftPrior");
        _muShiftPrior = settings->get<double>("muShiftPrior");
        _lambdaInitRootPrior = settings->get<double>("lambdaInitRootPrior");
        _muInitRootPrior = settings->get<double>("muInitRootPrior");
        _lambdaShiftRootPrior = settings->get<double>("lambdaShiftRootPrior");
        _muShiftRootPrior = settings->get<double>("muShiftRootPrior");
        _updateRateLambda0 = settings->get<double>("updateRateLambda0");
        _updateRateMu0 = settings->get<double>("updateRateMu0");
        _updateRateLambdaShift = settings->get<double>("updateRateLambdaShift");
        _updateRateMuShift = settings->get<double>("updateRateMuShift");
    } else if (modelType == "trait") {
        _betaInit = settings->get<double>("betaInit");
        _betaShiftInit = settings->get<double>("betaShiftInit");
        _betaInitPrior = settings->get<double>("betaInitPrior");
        _betaShiftPrior = settings->get<double>("betaShiftPrior");
        _betaInitRootPrior = settings->get<double>("betaInitRootPrior");
        _betaShiftRootPrior = settings->get<double>("betaShiftRootPrior");
        _updateRateBeta0 = settings->get<double>("updateRateBeta0");
        _updateRateBetaShift = settings->get<double>("updateRateBetaShift");
    }

    _poissonRatePrior = settings->get<double>("poissonRatePrior");
}


Prior::~Prior()
{
}


double Prior::lambdaShiftPrior(double x)
{
    return _rng->lnNormalPdf(0.0, _lambdaShiftPrior, x);
}


double Prior::generateLambdaShiftFromPrior()
{
    if (_updateRateLambdaShift <= _UPDATE_TOL) {
        return _lambdaShift0;
    } else {
        return _rng->normalRv(0.0, _lambdaShiftPrior);
    }
}


double Prior::lambdaInitPrior(double x)
{
    return _rng->lnExponentialPdf(_lambdaInitPrior, x);
}


double Prior::generateLambdaInitFromPrior()
{
    if (_updateRateLambda0 <= _UPDATE_TOL) {
        return _lambdaInit0;
    } else {
        return _rng->exponentialRv(_lambdaInitPrior);
    }
}


double Prior::muInitPrior(double x)
{
    return _rng->lnExponentialPdf(_muInitPrior, x);
}


double Prior::generateMuInitFromPrior()
{
    if (_updateRateMu0 <= _UPDATE_TOL) {
        return _muInit0;
    } else {
        return _rng->exponentialRv(_muInitPrior);
    }
}


double Prior::muShiftPrior(double x)
{
    return _rng->lnNormalPdf(0.0, _muShiftPrior, x);
}


double Prior::generateMuShiftFromPrior()
{
    if (_updateRateMuShift <= _UPDATE_TOL) {
        return _muShift0;
    } else {
        return _rng->normalRv(0.0, _muShiftPrior);
    }
}


double Prior::poissonRatePrior(double x)
{
    return _rng->lnExponentialPdf(_poissonRatePrior, x);
}


double Prior::generatePoissonRateFromPrior()
{
    return _rng->exponentialRv(_poissonRatePrior);
}


double Prior::betaInitPrior(double x)
{
    return _rng->lnExponentialPdf(_betaInitPrior, x);
}


double Prior::generateBetaInitFromPrior()
{
    if (_updateRateBeta0 <= _UPDATE_TOL) {
        return _betaInit;
    } else {
        return _rng->exponentialRv(_betaInitPrior);
    }
}


double Prior::betaShiftPrior(double x)
{
    return _rng->lnNormalPdf(0.0, _betaShiftPrior, x);
}


double Prior::generateBetaShiftFromPrior()
{
    if (_updateRateBetaShift <= _UPDATE_TOL) {
        return _betaShiftInit;
    } else {
        return _rng->normalRv(0.0, _betaShiftPrior);
    }
}


double Prior::lambdaInitRootPrior(double x)
{
    if (fabs(_lambdaInitRootPrior + 1) < _UPDATE_TOL) {
        return lambdaInitPrior(x);
    } else {
        return _rng->lnExponentialPdf(_lambdaInitRootPrior, x);
    }

}


double Prior::lambdaShiftRootPrior(double x)
{
    if (fabs(_lambdaShiftRootPrior + 1) < _UPDATE_TOL) {
        return lambdaShiftPrior(x);
    } else {
        return _rng->lnExponentialPdf(_lambdaShiftRootPrior, x);
    }
}


double Prior::muInitRootPrior(double x)
{
    if (fabs(_muInitRootPrior + 1) < _UPDATE_TOL) {
        return muInitPrior(x);
    } else {
        return _rng->lnExponentialPdf(_muInitRootPrior, x);
    }
}


double Prior::muShiftRootPrior(double x)
{
    if (fabs(_muShiftRootPrior + 1) < _UPDATE_TOL) {
        return muShiftPrior(x);
    } else {
        return _rng->lnExponentialPdf(_muShiftRootPrior, x);
    }
}


double Prior::betaInitRootPrior(double x)
{
    if (fabs(_betaInitRootPrior + 1) < _UPDATE_TOL) {
        return betaInitPrior(x);
    } else {
        return _rng->lnExponentialPdf(_betaInitRootPrior, x);
    }
}


double Prior::betaShiftRootPrior(double x)
{
    if (fabs(_betaShiftRootPrior + 1) < _UPDATE_TOL) {
        return betaShiftPrior(x);
    } else {
        return _rng->lnExponentialPdf(_betaShiftRootPrior, x);
    }
}
