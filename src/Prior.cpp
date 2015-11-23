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
        _lambdaIsTimeVariablePrior =
            settings->get<double>("lambdaIsTimeVariablePrior");
        _updateRateLambda0 = settings->get<double>("updateRateLambda0");
        _updateRateMu0 = settings->get<double>("updateRateMu0");
        _updateRateLambdaShift = settings->get<double>("updateRateLambdaShift");
        _updateRateMuShift = settings->get<double>("updateRateMuShift");
        
        /***************************/
        _preservationRatePrior = settings->get<double>("preservationRatePrior");
        
        
        
    } else if (modelType == "trait") {
        _betaInit = settings->get<double>("betaInit");
        _betaShiftInit = settings->get<double>("betaShiftInit");
        _betaInitPrior = settings->get<double>("betaInitPrior");
        _betaShiftPrior = settings->get<double>("betaShiftPrior");
        _betaInitRootPrior = settings->get<double>("betaInitRootPrior");
        _betaShiftRootPrior = settings->get<double>("betaShiftRootPrior");
        _betaIsTimeVariablePrior =
            settings->get<double>("betaIsTimeVariablePrior");
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
 
    
    // NOTE: bug discovered 10.25.2014, not sure when this was introduced.
    // was sending stdev^0.5 to lnNormalPDF rather than stdev
    //      (Stat::lnNormalPDF takes stdev as an argument).
    // The prior is a standard deviation, not variance
    
    // Bug version:
    //return Stat::lnNormalPDF(x, 0.0, std::sqrt(_lambdaShiftPrior));
    
    return Stat::lnNormalPDF(x, 0.0,  _lambdaShiftPrior);
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

double Prior::lambdaIsTimeVariablePrior()
{
    return _lambdaIsTimeVariablePrior;
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
 
    // Bug version: 10.25.2014
    // see above, lambdaShiftPrior
    //return Stat::lnNormalPDF(x, 0.0, std::sqrt(_muShiftPrior));
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


bool Prior::generateLambdaIsTimeVariableFromPrior()
{
    return _random.trueWithProbability(_lambdaIsTimeVariablePrior);
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

    
    // Bug version, fixed 10.25.2014
    //return Stat::lnNormalPDF(x, 0.0, std::sqrt(_betaShiftPrior));
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


bool Prior::generateBetaIsTimeVariableFromPrior()
{
    return _random.trueWithProbability(_betaIsTimeVariablePrior);
}

double Prior::betaIsTimeVariablePrior()
{
    return _betaIsTimeVariablePrior;
}


double Prior::lambdaInitRootPrior(double x)
{
    if (fabs(_lambdaInitRootPrior + 1) < _UPDATE_TOL) {
        return lambdaInitPrior(x);
    } else {
        return Stat::lnExponentialPDF(x, _lambdaInitRootPrior);
        std::cout << "Invalid call to lambdaInitRootPrior"  << std::endl;
        exit(0);
    }

}


double Prior::lambdaShiftRootPrior(double x)
{
    if (fabs(_lambdaShiftRootPrior + 1) < _UPDATE_TOL) {
        return lambdaShiftPrior(x);
    } else {
        //TODO: what is going on here?
        // lambdaShiftRootPrior has been deprecated,
        // but why is this implemented as an exponential???
        // TODO: throwing exception here as should never
        //  get here, but need more tests to remove.
        std::cout << "Invalid call to lambdaInitRootPrior"  << std::endl;
        exit(0);
        
        
        return Stat::lnExponentialPDF(x, _lambdaShiftRootPrior);
    }
}


double Prior::muInitRootPrior(double x)
{
    if (fabs(_muInitRootPrior + 1) < _UPDATE_TOL) {
        return muInitPrior(x);
    } else {

        // TODO: throwing exception here as should never
        //  get here, but need more tests to remove.
        std::cout << "Invalid call to muInitRootPrior"  << std::endl;
        exit(0);
        return Stat::lnExponentialPDF(x, _muInitRootPrior);
    }
}


double Prior::muShiftRootPrior(double x)
{
    if (fabs(_muShiftRootPrior + 1) < _UPDATE_TOL) {
        return muShiftPrior(x);
    } else {
        //TODO: what is going on here?
        // muShiftRootPrior has been deprecated,
        // but why is this implemented as an exponential???
        std::cout << "Invalid call to muShiftRootPrior"  << std::endl;
        exit(0);
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
        //TODO: what is going on here?
        // betaShiftRootPrior has been deprecated,
        // but why is this implemented as an exponential???
        return Stat::lnExponentialPDF(x, _betaShiftRootPrior);
    }
}


double Prior::preservationRatePrior(double x)
{

    //For now, fix a flat uniform prior on this.
    //Stat::lnExponentialPDF(x, _preservationRatePrior);
    return 0.0;
}

