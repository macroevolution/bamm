#include "MCMCDataWriter.h"
#include "Settings.h"
#include "Model.h"
#include "MCMC.h"

#include <iostream>


MCMCDataWriter::MCMCDataWriter(Settings& settings) :
    _outputFileName(settings.get("mcmcOutfile")),
    _outputFreq(settings.get<int>("mcmcWriteFreq"))
{
    initializeStream();
    writeHeader();
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
    return "generation,N_shifts,logPrior,logLik,eventRate,acceptRate";
}


MCMCDataWriter::~MCMCDataWriter()
{
    _outputStream.close();
}


void MCMCDataWriter::writeData(int generation, MCMC& mcmc)
{
    if (generation % _outputFreq != 0) {
        return;
    }

    Model& model = mcmc.model();

    _outputStream << generation                      << ","
                  << model.getNumberOfEvents()       << ","
                  << model.computeLogPrior()         << ","
                  << model.getCurrentLogLikelihood() << ","
                  << model.getEventRate()            << ","
                  << model.getMHAcceptanceRate()     << std::endl;
}
