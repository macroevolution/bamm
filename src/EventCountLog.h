//
//  EventCountLog.h
//  bammx
//
//  Created by Dan Rabosky on 5/23/14.
//  Copyright (c) 2014 Dan Rabosky. All rights reserved.
//

#ifndef __bammx__EventCountLog__
#define __bammx__EventCountLog__

#include <iostream>

class EventCountLog
{
    
public:
    
    EventCountLog(int x);
    
    
    int getAddProposeCount();
    void incrementAddProposeCount();
    int getAddAcceptCount();
    void incrementAddAcceptCount();
    int getSubtractProposeCount();
    void incrementSubtractProposeCount();
    int getSubtractAcceptCount();
    void incrementSubtractAcceptCount();
    
    int getEventCount();
 
    int getInStateCount();
    void incrementInStateCount();
    
private:

    int _addProposeCount;
    int _addAcceptCount;
    int _subtractProposeCount;
    int _subtractAcceptCount;
    
    int _eventCount;

    int _inStateCount;

};
 

inline int EventCountLog::getAddAcceptCount()
{
    return _addAcceptCount;
}

inline void EventCountLog::incrementAddAcceptCount()
{
    _addAcceptCount++;
}

inline int EventCountLog::getAddProposeCount()
{
    return _addProposeCount;
}

inline void EventCountLog::incrementAddProposeCount()
{
    _addProposeCount++;
}

inline int EventCountLog::getSubtractProposeCount()
{
    return _subtractProposeCount;
}

inline void EventCountLog::incrementSubtractProposeCount()
{
    _subtractProposeCount++;
}

inline int EventCountLog::getSubtractAcceptCount()
{
    return _subtractAcceptCount;
}

inline void EventCountLog::incrementSubtractAcceptCount()
{
    _subtractAcceptCount++;
}

inline int EventCountLog::getEventCount()
{
    return _eventCount;
}


inline int EventCountLog::getInStateCount()
{
    return _inStateCount;
}

inline void EventCountLog::incrementInStateCount()
{
    _inStateCount++;
}

#endif /* defined(__bammx__EventCountLog__) */
