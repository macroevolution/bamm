//
//  FastSimulatePrior.cpp
//  bamm
//
//  Created by Dan Rabosky on 11/29/13.
//  Copyright (c) 2013 Dan Rabosky. All rights reserved.
//

#include "FastSimulatePrior.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "Random.h"
#include "Settings.h"
#include "Stat.h"


FastSimulatePrior::FastSimulatePrior(Random& random, Settings* sp) :
    _random(random), sttings(sp)
{
    _generations = 0;

    _eventRate = 1 / sttings->get<double>("poissonRatePrior");
    
    _updateEventRateScale = sttings->get<double>("updateEventRateScale");
    
    _poissonRatePrior = sttings->get<double>("poissonRatePrior");
    _numberEvents = 0;
    
    _outfileName = sttings->get("priorOutputFileName");
    

    // Open streams for writing
    _fspOutStream.open(_outfileName.c_str());
 
    writeHeaderToOutputFile();
    
    std::cout << "\nSimulating prior distribution on shifts....\n" << std::endl;
    std::cout << "Progress: 0% ";
    std::cout.flush();
    
    int fints = (int)round(sp->get<int>("numberGenerations") / 32);
    
    
    for (int i = 0; i < sp->get<int>("numberGenerations"); i++){
                
        updateState();
    
        int mcmcWriteFreq = sp->get<int>("mcmcWriteFreq");
        if (mcmcWriteFreq > 0 && (i % mcmcWriteFreq) == 0){
            writeStateTofile();
        }
        
        if (fints > 0 && (i % fints) == 0){
            
            if (i == (fints*16)){
                std::cout << " 50% ";
            }
            
            std::cout << "|";
            std::cout.flush();
        }
        
    }
    
    std::cout << " 100% \n\nDone...Results written to file << " << _outfileName.c_str();
    std::cout << " >>\n" << std::endl;   
    
}

FastSimulatePrior::~FastSimulatePrior()
{
    _fspOutStream.close();
}

void FastSimulatePrior::writeStateTofile()
{
    writeStateToStream(_fspOutStream);
}

void FastSimulatePrior::writeStateToStream(std::ostream &outStream)
{
    outStream   << _generations      << ","
                << _numberEvents    << ","
                << _eventRate       << std::endl;
}

void FastSimulatePrior::updateEventRateMH()
{

    double oldEventRate = getEventRate();
    
    double cterm = exp( _updateEventRateScale * (_random.uniform() - 0.5) );
    setEventRate(cterm * oldEventRate);
    
    double LogPriorRatio =
        Stat::lnExponentialPDF(getEventRate(), _poissonRatePrior) -
        Stat::lnExponentialPDF(oldEventRate, _poissonRatePrior);
    double logProposalRatio = log(cterm);
    
    double logHR = LogPriorRatio + logProposalRatio;
    const bool acceptMove = acceptMetropolisHastings(logHR);
    
    
    if (acceptMove == false){
        setEventRate(oldEventRate);
    }

}


void FastSimulatePrior::changeNumberOfEventsMH()
{

    bool acceptMove = false;
    
    // Propose gains & losses equally if not on boundary (n = 0) events:
    
    // Current number of events on the tree, not counting root state:
    double K = (double)(_numberEvents);
    
    bool gain = _random.trueWithProbability(0.5);
    if (K == 0) {
        // set event to gain IF on boundary
        gain = true;
    }
    
    
    
    // now to adjust acceptance ratios:
    
    if (gain) {
        
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
    return _random.trueWithProbability(r);
}


void FastSimulatePrior::writeHeaderToOutputFile()
{
    _fspOutStream << "generation,N_shifts,eventRate" << std::endl;
}
