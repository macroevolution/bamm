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
 
    I have not implemented any of the generateFromPrior fxns
*/
    
    double lambdaShiftPrior(double);
    double generateLambdaShiftFromPrior(); // not done
    
    double lambdaInitPrior(double);
    double generateLambdaInitFromPrior(); // not done
    
    double muShiftPrior(double);
    double generateMuShiftFromPrior(); // not done
    
    double muInitPrior(double);
    double generateMuInitFromPrior(); // not done
    
    double poissonRatePrior(double);
    double generatePoissonRateFromPrior();
    
    double betaInitPrior(double);
    double generateBetaInitFromPrior();
    
    double betaShiftPrior(double);
    double generateBetaShiftFromPrior();
    
private:
    
    MbRandom *ranPtr;
    Settings *sttings;
};


#endif
