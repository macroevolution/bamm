#include "DataWriter.h"
#include "Settings.h"
#include "Model.h"

#include <iostream>
#include <iomanip>


DataWriter::DataWriter(Settings& settings, Model& model) :
    _model(model),
    _streamsInitialized(false),
    _headersWritten(false),
    _mcmcOutputFileName(settings.get("mcmcOutfile")), // TODO: Rename method
    _eventOutputFileName(settings.get("eventDataOutfile")),
    _acceptanceOutputFileName(settings.get("acceptanceInfoFileName")),
    _outputAcceptanceInfo(settings.get<bool>("outputAcceptanceInfo")),
    _mcmcOutputFreq(settings.get<int>("mcmcWriteFreq")),
    _eventOutputFreq(settings.get<int>("eventDataWriteFreq")),
    _stdoutOutputFreq(settings.get<int>("printFreq"))
{
}


DataWriter::~DataWriter()
{
    if (_streamsInitialized) {
        _mcmcOutputStream.close();
        _eventOutputStream.close();
        if (_outputAcceptanceInfo) {
            _acceptanceOutputStream.close();
        }
    }
}


void DataWriter::writeData(int generation)
{
    if (!_streamsInitialized) {
        initializeStreams();
        _streamsInitialized = true;
    }

    if (!_headersWritten) {
        writeHeaders();
        _headersWritten = true;
    }

    if (generation % _mcmcOutputFreq == 0) {
        writeMCMCData(generation);
    }

    if (generation % _eventOutputFreq == 0) {
        writeEventData(generation);
    }

    if (generation % _stdoutOutputFreq == 0) {
        writeStdoutData(generation);

        // Reset acceptance rates when state data are printed to the screen.
        _model.resetMHAcceptanceParameters();
    }

    writeAcceptanceData();
}


void DataWriter::initializeStreams()
{
    _mcmcOutputStream.open(_mcmcOutputFileName.c_str());
    _eventOutputStream.open(_eventOutputFileName.c_str());
    if (_outputAcceptanceInfo) {
        _acceptanceOutputStream.open(_acceptanceOutputFileName.c_str());
    }
}


void DataWriter::writeHeaders()
{
    writeMCMCHeaders();
    writeEventHeaders();
    writeStdoutHeaders();
    writeAcceptanceHeaders();
}


void DataWriter::writeMCMCHeaders()
{
    _mcmcOutputStream
        << "generation,N_shifts,logPrior,logLik,eventRate,acceptRate\n";
}


void DataWriter::writeEventHeaders()
{
    _eventOutputStream << "generation,leftchild,rightchild,abstime";
    writeSpecificEventHeaders();
}


void DataWriter::writeStdoutHeaders()
{
    std::cout << std::setw(15) << "Generation"
              << std::setw(15) << "LogLikelihood"
              << std::setw(15) << "NumShifts"
              << std::setw(15) << "LogPrior"
              << std::setw(15) << "AcceptRate\n";
}


void DataWriter::writeAcceptanceHeaders()
{
    _acceptanceOutputStream << "proposedParam,accepted\n";
}


void DataWriter::writeMCMCData(int generation)
{
    _mcmcOutputStream << generation                        << ","
                      << _model.getNumberOfEvents()       << ","
                      << _model.computeLogPrior()         << ","
                      << _model.getCurrentLogLikelihood() << ","
                      << _model.getEventRate()            << ","
                      << _model.getMHAcceptanceRate()     << std::endl;
}


void DataWriter::writeEventData(int generation)
{
    // TODO: Model shouldn't be concerned with output
    std::stringstream eventData;
    _model.getEventDataString(eventData, generation);
    _eventOutputStream << eventData.str() << std::endl;
}


void DataWriter::writeStdoutData(int generation)
{
    std::cout << std::setw(15) << generation
              << std::setw(15) << _model.getCurrentLogLikelihood()
              << std::setw(15) << _model.getNumberOfEvents()
              << std::setw(15) << _model.computeLogPrior()
              << std::setw(15) << _model.getMHAcceptanceRate() << std::endl;
}


void DataWriter::writeAcceptanceData()
{
    int param = _model.getLastParameterUpdated();
    int accepted = _model.getAcceptLastUpdate();
    _model.setAcceptLastUpdate(-1);

    _acceptanceOutputStream << param << "," << accepted << std::endl;
}
