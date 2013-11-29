//
//  FastSimulatePrior.h
//  bamm
//
//  Created by Dan Rabosky on 11/29/13.
//  Copyright (c) 2013 Dan Rabosky. All rights reserved.
//

#ifndef __bamm__FastSimulatePrior__
#define __bamm__FastSimulatePrior__

#include <iostream>

#include <vector>

//Forward declarations
class Node;
class MbRandom;
class Settings;

class FastSimulatePrior
{
public:
    
    FastSimulatePrior(MbRandom* ranptr, Settings* sp);
    ~FastSimulatePrior();
    
    void updateState(void);
    void changeNumberOfEventsMH(void);
    void updateEventRateMH(void);
    int getNumberOfEvents(void);
    double getEventRate(void);
    bool acceptMetropolisHastings(const double lnR);
    
    
private:

    MbRandom* ran;
    Settings* sttings;
    
    double _eventRate;
    int _numberEvents;
    int _generations;
    double _updateEventRateScale;
    double _poissonRatePrior;
    
    void incrementGeneration(void);
    void setEventRate(double x);
    

};

inline void FastSimulatePrior::setEventRate(double x)
{
    _eventRate = x;
}

inline void FastSimulatePrior::incrementGeneration()
{
    _generations++;
}

inline int FastSimulatePrior::getNumberOfEvents()
{
    return _numberEvents;
}


inline double FastSimulatePrior::getEventRate()
{
    return _eventRate;
}



#endif /* defined(__bamm__FastSimulatePrior__) */
