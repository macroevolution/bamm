#include "Prior.h"
#include "Random.h"
#include "Settings.h"
#include "Stat.h"

#define _UPDATE_TOL 0.0001


Prior::Prior(Random& random, Settings* settings) : _random(random)
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
    return Stat::lnNormalPDF(x, 0.0, _lambdaShiftPrior);
}


double Prior::generateLambdaShiftFromPrior()
{
    if (_updateRateLambdaShift <= _UPDATE_TOL) {
        return _lambdaShift0;
    } else {
        return _random.normal(0.0, _lambdaShiftPrior);
    }
}


double Prior::lambdaInitPrior(double x)
{
    return Stat::lnExponentialPDF(x, _lambdaInitPrior);
}


double Prior::generateLambdaInitFromPrior()
{
    if (_updateRateLambda0 <= _UPDATE_TOL) {
        return _lambdaInit0;
    } else {
        return _random.exponential(_lambdaInitPrior);
    }
}


double Prior::muInitPrior(double x)
{
    return Stat::lnExponentialPDF(x, _muInitPrior);
}


double Prior::generateMuInitFromPrior()
{
    if (_updateRateMu0 <= _UPDATE_TOL) {
        return _muInit0;
    } else {
        return _random.exponential(_muInitPrior);
    }
}


double Prior::muShiftPrior(double x)
{
    return Stat::lnNormalPDF(x, 0.0, _muShiftPrior);
}


double Prior::generateMuShiftFromPrior()
{
    if (_updateRateMuShift <= _UPDATE_TOL) {
        return _muShift0;
    } else {
        return _random.normal(0.0, _muShiftPrior);
    }
}


double Prior::poissonRatePrior(double x)
{
    return Stat::lnExponentialPDF(x, _poissonRatePrior);
}


double Prior::generatePoissonRateFromPrior()
{
    return _random.exponential(_poissonRatePrior);
}


double Prior::betaInitPrior(double x)
{
    return Stat::lnExponentialPDF(x, _betaInitPrior);
}


double Prior::generateBetaInitFromPrior()
{
    if (_updateRateBeta0 <= _UPDATE_TOL) {
        return _betaInit;
    } else {
        return _random.exponential(_betaInitPrior);
    }
}


double Prior::betaShiftPrior(double x)
{
    return Stat::lnNormalPDF(x, 0.0, _betaShiftPrior);
}


double Prior::generateBetaShiftFromPrior()
{
    if (_updateRateBetaShift <= _UPDATE_TOL) {
        return _betaShiftInit;
    } else {
        return _random.normal(0.0, _betaShiftPrior);
    }
}


double Prior::lambdaInitRootPrior(double x)
{
    if (fabs(_lambdaInitRootPrior + 1) < _UPDATE_TOL) {
        return lambdaInitPrior(x);
    } else {
        return Stat::lnExponentialPDF(x, _lambdaInitRootPrior);
    }

}


double Prior::lambdaShiftRootPrior(double x)
{
    if (fabs(_lambdaShiftRootPrior + 1) < _UPDATE_TOL) {
        return lambdaShiftPrior(x);
    } else {
        return Stat::lnExponentialPDF(x, _lambdaShiftRootPrior);
    }
}


double Prior::muInitRootPrior(double x)
{
    if (fabs(_muInitRootPrior + 1) < _UPDATE_TOL) {
        return muInitPrior(x);
    } else {
        return Stat::lnExponentialPDF(x, _muInitRootPrior);
    }
}


double Prior::muShiftRootPrior(double x)
{
    if (fabs(_muShiftRootPrior + 1) < _UPDATE_TOL) {
        return muShiftPrior(x);
    } else {
        return Stat::lnExponentialPDF(x, _muShiftRootPrior);
    }
}


double Prior::betaInitRootPrior(double x)
{
    if (fabs(_betaInitRootPrior + 1) < _UPDATE_TOL) {
        return betaInitPrior(x);
    } else {
        return Stat::lnExponentialPDF(x, _betaInitRootPrior);
    }
}


double Prior::betaShiftRootPrior(double x)
{
    if (fabs(_betaShiftRootPrior + 1) < _UPDATE_TOL) {
        return betaShiftPrior(x);
    } else {
        return Stat::lnExponentialPDF(x, _betaShiftRootPrior);
    }
}
