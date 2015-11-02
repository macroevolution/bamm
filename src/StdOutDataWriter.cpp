#include "StdOutDataWriter.h"
#include "Settings.h"
#include "Model.h"

#include <iostream>
#include <iomanip>


StdOutDataWriter::StdOutDataWriter(Settings& settings) :
    _outputFreq(settings.get<int>("printFreq")),
    _headerWritten(false)
{
}


StdOutDataWriter::~StdOutDataWriter()
{
}


void StdOutDataWriter::writeData(int generation, Model& model)
{
    if (!_headerWritten && _outputFreq > 0) {
        writeHeader();
        _headerWritten = true;
    }

    if (_outputFreq == 0 || generation % _outputFreq != 0) {
        return;
    }

    std::cout << std::setw(12) << generation
              << std::setw(12) << model.getNumberOfEvents()
              << std::setw(12) << model.computeLogPrior()
              << std::setw(12) << model.getCurrentLogLikelihood()
              << std::setw(12) << model.getEventRate()
              << std::setw(12) << model.getMHAcceptanceRate()
              << std::endl;
}


void StdOutDataWriter::writeHeader()
{
    std::cout << header() << std::endl;
}


std::string StdOutDataWriter::header()
{
    return "  generation"
           "    N_shifts"
           "    logPrior"
           "      logLik"
           "   eventRate"
           "  acceptRate";
}
