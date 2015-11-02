#ifndef CHAIN_SWAP_DATA_WRITER_H
#define CHAIN_SWAP_DATA_WRITER_H


#include <vector>
#include <string>
#include <fstream>

class Settings;
class Model;
class MCMC;


class ChainSwapDataWriter
{
public:

    ChainSwapDataWriter(Settings& settings);
    ~ChainSwapDataWriter();

    void writeData(int generation, const std::vector<MCMC*>& chains,
        int chain_1, int chain_2, bool accepted);

private:

    void initializeStream();
    void writeHeader();
    std::string header() const;

    std::vector<int> rankChainsByTemp(const std::vector<MCMC*>& chains) const;
    std::vector<double> chainTemperatures
        (const std::vector<MCMC*>& chains) const;
    std::vector<double> sortValues(std::vector<double> values) const;
    int rankValue(double value, std::vector<double> sortedValues) const;

    int _numberOfChains;

    std::string _outputFileName;
    std::ofstream _outputStream;
};


#endif
