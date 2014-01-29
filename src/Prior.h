//
//  Prior.h
//  bamm
//
//  Created by Dan Rabosky on 1/28/14.
//  Copyright (c) 2014 Dan Rabosky. All rights reserved.
//

#ifndef __bamm__Prior__
#define __bamm__Prior__

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
    double generateLambdaShiftFromPrior(void); // not done
    
    double lambdaInitPrior(double);
    double generateLambdaInitFromPrior(void); // not done
    
    double muShiftPrior(double);
    double generateMuShiftFromPrior(void); // not done
    
    double muInitPrior(double);
    double generateMuInitFromPrior(double); // not done
    
    double poissonRatePrior(double);
    double generatePoissonRateFromPrior(double);
    
private:
    
    MbRandom *ranPtr;
    Settings *sttings;
    
    

};



#endif /* defined(__bamm__Prior__) */
