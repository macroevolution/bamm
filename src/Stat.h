#ifndef STAT_H
#define STAT_H


#include "MbRandom.h"
#include <vector>

#define MAX_LAG 1000

class Stat
{
public:

    static double standard_deviation(const std::vector<double>& values);
    static double variance(const std::vector<double>& values);

    static double lnNormalPDF(double x, double mean, double sd);
    static double lnExponentialPDF(double x, double rate);

    static double ESS(const std::vector<double>& values, double burninFrac);

private:

    static MbRandom _random;
};


#endif
