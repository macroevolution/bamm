//
//  EventCountLog.cpp
//  bammx
//
//  Created by Dan Rabosky on 5/23/14.
//  Copyright (c) 2014 Dan Rabosky. All rights reserved.
//

#include "EventCountLog.h"


EventCountLog::EventCountLog(int x)
{
    _eventCount = x;

    _addProposeCount = 0;
    _addAcceptCount = 0;
    _subtractProposeCount = 0;
    _subtractAcceptCount = 0;

    _inStateCount = 0;

}


