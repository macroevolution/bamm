#include "AcceptanceDataWriter.h"
#include "Settings.h"
#include "Model.h"
#include "MCMC.h"

#include <iostream>


AcceptanceDataWriter::AcceptanceDataWriter(Settings& settings) :
    _shouldOutputData(settings.get<bool>("outputAcceptanceInfo")),
    _outputFileName(settings.get("acceptanceInfoFileName"))
{
    if (_shouldOutputData) {
        initializeStream();
        writeHeader();
    }
}


void AcceptanceDataWriter::initializeStream()
{
    _outputStream.open(_outputFileName.c_str());
}


void AcceptanceDataWriter::writeHeader()
{
    _outputStream << header() << std::endl;
}


std::string AcceptanceDataWriter::header()
{
    return "param,accepted";
}


AcceptanceDataWriter::~AcceptanceDataWriter()
{
    if (_shouldOutputData) {
        _outputStream.close();
    }
}


void AcceptanceDataWriter::writeData(MCMC& mcmc)
{
    if (!_shouldOutputData) {
        return;
    }

    Model& model = mcmc.model();
    _outputStream << model.getLastParameterUpdated() << ","
                  << model.getAcceptLastUpdate()     << std::endl;
}
