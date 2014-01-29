//
//  Prior.cpp
//  bamm
//
//  Created by Dan Rabosky on 1/28/14.
//  Copyright (c) 2014 Dan Rabosky. All rights reserved.
//

#include "Prior.h"
#include "MbRandom.h"
#include "Settings.h"

Prior::Prior(MbRandom* ran, Settings* sp)
{
    
    ranPtr = ran;
    sttings = sp;
    
    std::cout << "Initialized class Prior object\n\n" << std::endl;

}

Prior::~Prior()
{

}

double Prior::lambdaShiftPrior(double x){

    double pp = 0;
    pp = ranPtr->lnNormalPdf((double)0.0, sttings->getLambdaShiftPrior(), x);
    
    
    return pp;
}


double Prior::lambdaInitPrior(double x){

    double pp = ranPtr->lnExponentialPdf(sttings->getLambdaInitPrior(), x);
    
    return pp;
    
}


double Prior::muInitPrior(double x){
    
    double pp = ranPtr->lnExponentialPdf(sttings->getMuInitPrior(), x);
    
    return pp;

}


double Prior::muShiftPrior(double x){

    double pp = ranPtr->lnNormalPdf((double)0.0, sttings->getMuShiftPrior(), x);
    
    return pp;

}



double Prior::poissonRatePrior(double x){
     
    double pp = ranPtr->lnExponentialPdf(sttings->getPoissonRatePrior(), x);
    
    return pp;

}















