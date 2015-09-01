//
//  FossilData.h
//  
//
//  Created by Dan Rabosky on 8/29/15.
//
//

#ifndef ____FossilData__
#define ____FossilData__


#include <stdio.h>
#include <string>
#include <vector>

class FossilData
{
    
public:
    FossilData(std::string fname);
    ~FossilData();
    
private:
    
    
    std::vector<std::string> _stagenames;
    std::vector<double> _startTime;
    std::vector<double> _endTime;
    std::vector<double> _relPresRate;
    std::vector<double> _fossilCount;    
    
    
};




#endif /* defined(____FossilData__) */




