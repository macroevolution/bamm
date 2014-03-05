#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "Autotune.h"

#include "SpExModel.h"
#include "TraitModel.h"
#include "MbRandom.h"
#include "Node.h"
#include "Tree.h"
#include "Settings.h"
#include "BranchEvent.h"
#include "TraitBranchEvent.h"

Autotune::Autotune(MbRandom* ran, SpExModel * mymodel, Settings*sp)
{
    std::cout << "\n\nThis is a tuning run to optimize MCMC operators\n" << std::endl;
    
    ranPtr = ran;
    ModelPtr = mymodel;
    sttings = sp;    
    std::string tunefilename = "autotune.csv"; // set this in settings instead
    
    std::ifstream outStream(tunefilename.c_str());
    
    bool fileOverwrite = sttings->getOverwrite();

    if (!fileOverwrite) {
        if (outStream) {
            exitWithErrorOutputFileExists(tunefilename);
        }
    }
    
    _lowerBound_EventRateScale = 0.0;
    _upperBound_EventRateScale = 80.0;
    _lowerBound_MoveScale = 0;
    _upperBound_MoveScale = 20;
    _lowerBound_lambdaInitScale = 0.0;
    _upperBound_lambdaInitScale = 10.0;
    _lowerBound_lambdaShiftScale = 0.0;
    _upperBound_lambdaShiftScale = 0.3;
    _lowerBound_muInitScale = 0.0;
    _upperBound_muInitScale = 10.0;
    _lowerBound_muShiftScale = 0.0;
    _upperBound_muShiftScale = 1.0;
    _lowerBound_betaInitScale = 0.0;
    _upperBound_betaInitScale = 10.0;
    _lowerBound_betaShiftScale = 0.0;
    _upperBound_betaShiftScale = 1.0;
    
    int TUNE_GENS = 200000;    
    int PRINT_FREQ = 1000;
    
    setUpdateWeights_Diversification();
    ModelPtr->resetGeneration();
    
    for (int i = 0; i < TUNE_GENS; i++) {
        int parmToUpdate = pickParameterClassToUpdate();
        updateState_Diversification(parmToUpdate); // update state        
        
        if ((i % PRINT_FREQ) == 0) {
            std::cout << std::setw(15) << "Tune_generation" << std::setw(15) << i << std::endl;
        }
    }
    
    std::ofstream ss;
    ss.open(tunefilename.c_str(), std::ofstream::app);
    for (int i = 0; i < TUNE_GENS; i++) { 
        ss << _param_ID[i] << "," << _param_value[i] << ",";
        ss << _accept_move[i] << "\n";
    }
    ss.close();
}


Autotune::~Autotune(void)
{

}



void Autotune::setUpdateWeights_Diversification(void)
{
    parWts.push_back(sttings->getUpdateRateEventNumber()); // event number
    parWts.push_back(sttings->getUpdateRateEventPosition()); // event position
    parWts.push_back(sttings->getUpdateRateEventRate()); // event rate
    parWts.push_back(sttings->getUpdateRateLambda0()); // lambda0 rate
    parWts.push_back(sttings->getUpdateRateLambdaShift()); // lambda shift
    parWts.push_back(sttings->getUpdateRateMu0()); //mu0 rate
    parWts.push_back(sttings->getUpdateRateMuShift()); // mu shift
    
    double sumwts = parWts[0];
    for (std::vector<double>::size_type i = 1; i < parWts.size(); i++) {
        sumwts += parWts[i];
        parWts[i] += parWts[i - 1];
    }
    
    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++) {
        parWts[i] /= sumwts;
    }
}


int Autotune::pickParameterClassToUpdate(void)
{
    double rn = ranPtr->uniformRv();
    int parm = 0;
    for (std::vector<double>::size_type i = 0; i < parWts.size(); i++) {
        if (rn <= parWts[i]) {
            parm = (int)i;
            break;
        }
    }
    return parm;
}


void Autotune::updateState_Diversification(int parm)
{
    _param_ID.push_back(parm);
    
    if (parm == 0) {
        _param_value.push_back(-1.0);
        ModelPtr->changeNumberOfEventsMH();    
    } else if (parm == 1) {
        
        double nv = ranPtr->uniformRv(_lowerBound_MoveScale, _upperBound_MoveScale);
        ModelPtr->setMoveSizeScale(nv);
        _param_value.push_back(nv);
        ModelPtr->moveEventMH();
        
    } else if (parm == 2) {
        
        double nv = ranPtr->uniformRv(_lowerBound_EventRateScale, _upperBound_EventRateScale);
        ModelPtr->setUpdateEventRateScale(nv);
        ModelPtr->updateEventRateMH();
        _param_value.push_back(nv);
        
    } else if (parm == 3) {
        double nv = ranPtr->uniformRv(_lowerBound_lambdaInitScale, _upperBound_lambdaInitScale);
        ModelPtr->setUpdateLambdaInitScale(nv);
        ModelPtr->updateLambdaInitMH();
        _param_value.push_back(nv);
        
    } else if (parm == 4) {
        double nv = ranPtr->uniformRv(_lowerBound_lambdaShiftScale, _upperBound_lambdaShiftScale);
        ModelPtr->setUpdateLambdaShiftScale(nv);
        ModelPtr->updateLambdaShiftMH();
        _param_value.push_back(nv);
        
    } else if (parm == 5) {
        
        double nv=ranPtr->uniformRv(_lowerBound_muInitScale, _upperBound_muInitScale);
        ModelPtr->setUpdateMuInitScale(nv);
        ModelPtr->updateMuInitMH();
        _param_value.push_back(nv);
        
    } else if (parm == 6) {
    
        std::cout << "Does not support autotuning of mu_shift param yet" << std::endl;
        throw;
        
    } else if (parm == 7) {
    
        std::cout << "Unsupported option (7) in" << std::endl;
        std::cout <<"autotune::updateState_Diversification" << std::endl;
        throw;
         
    } else {
        // should never get here...throw exception?
        std::cout << "Bad parm to update\n" << std::endl;
    }
    
    _accept_move.push_back(ModelPtr->getAcceptLastUpdate());
 
    // reset to unmodified value
    ModelPtr->setAcceptLastUpdate(-1);
}

void Autotune::exitWithErrorOutputFileExists(std::string const& tunefilename)
{
    std::cout << "ERROR: Analysis is set to not overwrite files.\n";
    std::cout << "Existing file '" << tunefilename << "' will be overwritten.\n";
    std::cout << "Fix by removing or renaming output file,\n";
    std::cout << "or set \"overwrite = 1\" in the control file.\n";
    std::exit(1);
}

