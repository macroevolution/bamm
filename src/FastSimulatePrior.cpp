//
//  FastSimulatePrior.cpp
//  bamm
//
//  Created by Dan Rabosky on 11/29/13.
//  Copyright (c) 2013 Dan Rabosky. All rights reserved.
//

#include "FastSimulatePrior.h"


#include <cmath>

#include "MbRandom.h"
#include "Settings.h"

FastSimulatePrior::FastSimulatePrior(MbRandom* ranptr, Settings* sp)
{
    
    ran = ranptr;    
    sttings = sp;
    
    for (int i =0; i <= 1000; i++){
        ran->uniformRv();
    }
    
    

    _eventRate = 1 / sttings->getPoissonRatePrior();
    
    _updateEventRateScale = sttings->getUpdateEventRateScale();
    
    _poissonRatePrior = sttings->getPoissonRatePrior();
    _numberEvents = 0;
    
}

FastSimulatePrior::~FastSimulatePrior()
{

}


void FastSimulatePrior::updateEventRateMH()
{

    double oldEventRate = getEventRate();
    
    double cterm = exp( _updateEventRateScale * (ran->uniformRv() - 0.5) );
    setEventRate(cterm * oldEventRate);
    
    double LogPriorRatio = ran->lnExponentialPdf(_poissonRatePrior,
                                                 getEventRate()) - ran->lnExponentialPdf(_poissonRatePrior, oldEventRate);
    double logProposalRatio = log(cterm);
    
    
    double logHR = LogPriorRatio + logProposalRatio;
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    
    if (acceptMove == false){
        setEventRate(oldEventRate);
    }

}


void FastSimulatePrior::changeNumberOfEventsMH()
{

    int proposedState = 0;
    bool acceptMove = false;
    
    // Propose gains & losses equally if not on boundary (n = 0) events:
    
    // Current number of events on the tree, not counting root state:
    double K = (double)(_numberEvents);
    
    bool gain = (ran->uniformRv() <= 0.5);
    if (K == 0) {
        // set event to gain IF on boundary
        gain = true;
    }
    
    
    
    // now to adjust acceptance ratios:
    
    if (gain) {
        
        proposedState = _numberEvents + 1;
        
        double qratio = 1.0;
        if (K == 0) {
            // no events on tree
            // can only propose gains.
            qratio = 0.5;
        } else {
            // DO nothing.
        }
        
        // Prior ratio is eventRate / (k + 1)
        // but now, eventCollection.size() == k + 1
        //  because event has already been added.
        // Here HR is just the prior ratio
        
        double logHR = log(_eventRate) - log(K + 1.0);
        
        // Now add log qratio
        
        logHR += log(qratio);
        acceptMove = acceptMetropolisHastings(logHR);
        
        if (acceptMove) {
            
            _numberEvents++;
            
        } 

    } else {

        double qratio = 1.0; // if loss, can only be qratio of 1.0
        if (K  == 1)
            qratio = 2.0;
        
        // This is probability of going from k to k-1
        // So, prior ratio is (k / eventRate)
        
        // First get prior ratio:
        double logHR = log(K) - log(_eventRate);
        
        // Now correct for proposal ratio:
        logHR += log(qratio);
        
        acceptMove = acceptMetropolisHastings(logHR);
        
        if (acceptMove) {
            
            _numberEvents--;

        }

    }
    
}


void FastSimulatePrior::updateState()
{

    updateEventRateMH();
    changeNumberOfEventsMH();
    
    incrementGeneration();

}


bool FastSimulatePrior::acceptMetropolisHastings(const double lnR)
{
    const double r = exp(lnR);
    return (ran->uniformRv() < r);
}





