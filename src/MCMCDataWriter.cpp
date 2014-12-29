#include "MCMCDataWriter.h"
#include "Settings.h"
#include "Model.h"

// TODO: make abstract class for MCMC datawriter with derived classes
//                              than can handle fossils
//  Then no more need to declare SpExModel.h here.

#include "SpExModel.h"

#include <iostream>


MCMCDataWriter::MCMCDataWriter(Settings& settings) :
    _outputFileName(settings.get("mcmcOutfile")),
    _outputFreq(settings.get<int>("mcmcWriteFreq"))
{
    
    if (settings.get<double>("updateRatePreservationRate") >= 0.0){
        _hasPreservationRate = true;
    }else{
        _hasPreservationRate = false;
    }
    
    
    if (_outputFreq > 0) {
        initializeStream();
        writeHeader();
    }
    

}


void MCMCDataWriter::initializeStream()
{
    _outputStream.open(_outputFileName.c_str());
}


void MCMCDataWriter::writeHeader()
{
    _outputStream << header() << std::endl;
}


std::string MCMCDataWriter::header()
{
    if (_hasPreservationRate){
        return "generation,N_shifts,logPrior,logLik,eventRate,preservationRate,acceptRate";
    }else{
        return "generation,N_shifts,logPrior,logLik,eventRate,acceptRate";    
    }
}


MCMCDataWriter::~MCMCDataWriter()
{
    if (_outputFreq > 0) {
        _outputStream.close();
    }
}


void MCMCDataWriter::writeData(int generation, Model& model)
{
    if (_outputFreq == 0 || generation % _outputFreq != 0) {
        return;
    }

    if (_hasPreservationRate == false){
        _outputStream << generation                        << ","
                  << model.getNumberOfEvents()             << ","
                  << model.computeLogPrior()               << ","
                  << model.getCurrentLogLikelihood()       << ","
                  << model.getEventRate()                  << ","
                  << model.getMHAcceptanceRate()           << std::endl;
    }else{
       
        double prate = static_cast<SpExModel*>(&model)->getPreservationRate();
        
        _outputStream << generation                        << ","
        << model.getNumberOfEvents()                       << ","
        << model.computeLogPrior()                         << ","
        << model.getCurrentLogLikelihood()                 << ","
        << model.getEventRate()                            << ","
        << prate    << ","
        << model.getMHAcceptanceRate()                     << std::endl;
    
    }
    

}




