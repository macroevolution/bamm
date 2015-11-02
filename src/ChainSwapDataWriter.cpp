#include "ChainSwapDataWriter.h"
#include "Settings.h"
#include "Model.h"
#include "MCMC.h"

#include <iostream>
#include <vector>
#include <algorithm>


ChainSwapDataWriter::ChainSwapDataWriter(Settings& settings) :
    _numberOfChains(settings.get<int>("numberOfChains")),
    _outputFileName(settings.get("chainSwapFileName"))
{
    if (_numberOfChains > 1) {
        initializeStream();
        writeHeader();
    }
}


void ChainSwapDataWriter::initializeStream()
{
    _outputStream.open(_outputFileName.c_str());
}


void ChainSwapDataWriter::writeHeader()
{
    _outputStream << header() << std::endl;
}


std::string ChainSwapDataWriter::header() const
{
    return "generation,rank_1,rank_2,swapAccepted";
}


ChainSwapDataWriter::~ChainSwapDataWriter()
{
    if (_numberOfChains > 1) {
        _outputStream.close();
    }
}


void ChainSwapDataWriter::writeData(int generation,
    const std::vector<MCMC*>& chains, int chain_1, int chain_2, bool accepted)
{
    const std::vector<int>& chainRanks = rankChainsByTemp(chains);

    int rank_1 = chainRanks[chain_1];
    int rank_2 = chainRanks[chain_2];

    if (rank_2 < rank_1) {
        std::swap(rank_1, rank_2);
    }

    _outputStream << generation  << ","
                  << rank_1      << ","
                  << rank_2      << ","
                  << accepted    << std::endl;
}


std::vector<int> ChainSwapDataWriter::rankChainsByTemp
    (const std::vector<MCMC*>& chains) const
{
    const std::vector<double>& temps = chainTemperatures(chains);
    const std::vector<double>& sortedTemps = sortValues(temps);

    std::vector<int> ranks;
    for (int i = 0; i < (int)temps.size(); i++) {
        int tempRank = rankValue(temps[i], sortedTemps);
        if (tempRank > 0) {
            ranks.push_back(tempRank);
        } else {
            log(Error) << "Error while ranking chain temperatures.\n";
            std::exit(1);
        }
    }

    return ranks;
}


std::vector<double> ChainSwapDataWriter::chainTemperatures
    (const std::vector<MCMC*>& chains) const
{
    std::vector<double> temps;
    for (int i = 0; i < (int)chains.size(); i++) {
        temps.push_back(chains[i]->model().getTemperatureMH());
    }
    return temps;
}


std::vector<double> ChainSwapDataWriter::sortValues
    (std::vector<double> values) const
{
    std::vector<double> sortedValues(values);
    std::sort(sortedValues.begin(), sortedValues.end());
    std::reverse(sortedValues.begin(), sortedValues.end());
    return sortedValues;
}


int ChainSwapDataWriter::rankValue
    (double value, std::vector<double> sortedValues) const
{
    for (int i = 0; i < (int)sortedValues.size(); i++) {
        if (value == sortedValues[i])
            return i + 1;
    }

    return -1;
}
