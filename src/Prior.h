//  Returns either:
//      Log-prior densities for parameter
//      A random deviate from the prior distribution for a parameter

#ifndef PRIOR_H
#define PRIOR_H

#include <iostream>

class Settings;
class Random;


class Prior
{
public:
    
    Prior(Random& random, Settings* sp);
    ~Prior();

/*
    Every function to compute a prior density should have a 
    counterpart function that generates a random variable from
    the same distribution.
 
*/
    
    double lambdaShiftPrior(double);
    double generateLambdaShiftFromPrior();
    
    double lambdaInitPrior(double);
    double generateLambdaInitFromPrior(); 
    
    double muShiftPrior(double);
    double generateMuShiftFromPrior(); 
    
    double muInitPrior(double);
    double generateMuInitFromPrior(); 

    bool generateLambdaIsTimeVariableFromPrior();
    double lambdaIsTimeVariablePrior();
    
    double poissonRatePrior(double);
    double generatePoissonRateFromPrior();
    
    double betaInitPrior(double);
    double generateBetaInitFromPrior();
    
    double betaShiftPrior(double);
    double generateBetaShiftFromPrior();

    bool generateBetaIsTimeVariableFromPrior();
    double betaIsTimeVariablePrior();
    
    double preservationRatePrior(double);
    
    
// Root priors:
    
    double lambdaShiftRootPrior(double);
    double generateLambdaShiftRootPrior();

    double lambdaInitRootPrior(double);
    double generateLambdaInitRootFromPrior();
    
    double muShiftRootPrior(double);
    double generateMuShiftRootFromPrior();
    
    double muInitRootPrior(double);
    double generateMuInitRootFromPrior();
    
    double betaInitRootPrior(double);
    double generateBetaInitRootFromPrior();
    
    double betaShiftRootPrior(double);
    double generateBetaShiftRootFromPrior();
    

private:
    
    Random& _random;

    // Initial parameters
    double _lambdaInit0;
    double _muInit0;
    double _betaInit;

    // Shift parameters
    double _lambdaShift0;
    double _muShift0;
    double _betaShiftInit;

    // Priors (non-root) for initial parameters
    double _lambdaInitPrior;
    double _muInitPrior;
    double _betaInitPrior;

    // Priors (non-root) for shift parameters
    double _lambdaShiftPrior;
    double _muShiftPrior;
    double _betaShiftPrior;

    // Root priors for initial parameters
    double _lambdaInitRootPrior;
    double _muInitRootPrior;
    double _betaInitRootPrior;

    // Root priors for shift parameters
    double _lambdaShiftRootPrior;
    double _muShiftRootPrior;
    double _betaShiftRootPrior;

    // Time variable/constant prior
    double _lambdaIsTimeVariablePrior;
    double _betaIsTimeVariablePrior;

    // Update rate for initial parameters
    double _updateRateLambda0;
    double _updateRateMu0;
    double _updateRateBeta0;

    // Update rate for shift parameters
    double _updateRateLambdaShift;
    double _updateRateMuShift;
    double _updateRateBetaShift;

    double _poissonRatePrior;
    
    double _preservationRatePrior;
    
};


#endif
