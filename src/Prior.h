//  Returns either:
//      Log-prior densities for parameter
//      A random deviate from the prior distribution for a parameter

#ifndef PRIOR_H
#define PRIOR_H

#include <iostream>

class Settings;
class MbRandom;


class Prior
{
public:
    
    Prior(MbRandom* ran, Settings* sp);
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
    
    double poissonRatePrior(double);
    double generatePoissonRateFromPrior();
    
    double betaInitPrior(double);
    double generateBetaInitFromPrior();
    
    double betaShiftPrior(double);
    double generateBetaShiftFromPrior();
    
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
    
    MbRandom *ranPtr;
    Settings *sttings;
};


#endif
