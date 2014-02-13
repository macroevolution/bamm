 

#ifndef Autotune_H 
#define Autotune_H
 

#include <stdlib.h>
#include <string>
#include <vector>
#include <iosfwd>

class MbRandom;
class SpExModel;
class TraitModel;
class Settings;


class Autotune
{
public:    
    
    Autotune(MbRandom* ran, SpExModel* mymodel, Settings* sp);
    ~Autotune();
    
    void setUpdateWeights_Diversification(void);
    int pickParameterClassToUpdate(void);
    void updateState_Diversification(int parm);
    void exitWithErrorOutputFileExists(std::string const& tunefilename);
    
private:    
    MbRandom* ranPtr;
    SpExModel*    ModelPtr;
    Settings* sttings;
    
    std::vector<double> parWts;    
    
    std::vector<int> _param_ID;
    std::vector<int> _accept_move;
    std::vector<double> _param_value;
    
// bounds for tuning:
    double _lowerBound_EventRateScale;
    double _upperBound_EventRateScale;
    double _lowerBound_MoveScale;
    double _upperBound_MoveScale;
    double _lowerBound_lambdaInitScale;
    double _upperBound_lambdaInitScale;
    double _lowerBound_lambdaShiftScale;
    double _upperBound_lambdaShiftScale;
    double _lowerBound_muInitScale;
    double _upperBound_muInitScale;
    double _lowerBound_muShiftScale;
    double _upperBound_muShiftScale;
    double _lowerBound_betaInitScale;
    double _upperBound_betaInitScale;
    double _lowerBound_betaShiftScale;
    double _upperBound_betaShiftScale;
    
};


#endif
