//
//  FastSimulatePrior.h
//  bamm
//
//  Created by Dan Rabosky on 11/29/13.
//  Copyright (c) 2013 Dan Rabosky. All rights reserved.
//

#ifndef __bamm__FastSimulatePrior__
#define __bamm__FastSimulatePrior__

#include <string>
#include <fstream>
#include <iosfwd>
#include <vector>
#include <cmath>

//Forward declarations
class Node;
class Random;
class Settings;
class EventCountLog;

class FastSimulatePrior
{
public:
    
    FastSimulatePrior(Random& random, Settings* sp);
    ~FastSimulatePrior();
    
    void updateState();
    void updateState(int min, int max);
    void changeNumberOfEventsMH();
    void changeNumberOfEventsMH(int min, int max);
    
    void updateEventRateMH();
    int getNumberOfEvents();
    double getEventRate();
    bool acceptMetropolisHastings(const double lnR);
    

    void fastSimulatePriorOldWay();
    void fastSimulatePriorExperimental();
    
// Output settings:
//    void writeStateToFile();
    void writeHeaderToOutputFile();
    
private:

    int round(double x);

    Random& _random;
    Settings* sttings;
    
    double _eventRate;
    int _numberEvents;
    int _generations;
    double _updateEventRateScale;
    double _poissonRatePrior;
    
    std::string _outfileName;
    
    void incrementGeneration();
    void setEventRate(double x);

    std::ofstream _fspOutStream;
    std::ofstream _eventDataOutStream;
    
    void writeStateTofile();
    void writeStateToStream(std::ostream& outStream);
    void exitWithErrorOutputFileExists();

    
    
    /**** New params May 23 2014 ****/
    
    // Track event proposals: additions and subtractions.
    std::vector<EventCountLog*> _TrackingVector;
    
    void writeTransitionProbsToFile();
    void writePriorProbsToFile_Experimental();
    void writePriorProbsToFile_OldWay();

    int _maxEvents;
    int _intervalGens;
    
    // End new params
    
    
};


inline int FastSimulatePrior::round(double x)
{
    return std::ceil(x - 0.5);
}


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
