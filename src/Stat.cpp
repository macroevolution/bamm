#include "Stat.h"
#include "MbRandom.h"

#include <vector>
#include <cmath>
#include <cstddef>


MbRandom Stat::_random;


double Stat::standard_deviation(const std::vector<double>& values)
{
    return std::sqrt(variance(values));
}


// Algorithm taken from the Naive algorithm in Wikipedia:
// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
double Stat::variance(const std::vector<double>& values)
{
    int n = 0;
    double sum = 0.0;
    double sum_of_squares = 0.0;

    std::vector<double>::const_iterator it;
    for (it = values.begin(); it != values.end(); ++it) {
        n++;
        sum += *it;
        sum_of_squares += *it * *it;
    }

    return (sum_of_squares - sum * sum / (double)n) / (double)(n - 1);
}


double Stat::lnNormalPDF(double x, double mean, double sd)
{
    return _random.lnNormalPdf(mean, sd * sd, x);
}


double Stat::lnExponentialPDF(double x, double rate)
{
    return _random.lnExponentialPdf(rate, x);
}


// ESS code modified from RevBayes TraceAnalysisContinuous.cpp
// originally written by Sebastian Hoehna, GPL License
double Stat::ESS(const std::vector<double>& values, double burninFrac)
{
    if (values.size() <= 0) {
        return 0;
    }


    double m = 0;
    size_t size = values.size();
    size_t burnin = floor(values.size() * burninFrac);
    for (size_t ii = burnin; ii < size; ii++) {
        m += values.at(ii);
    }

    double mean = m / double(size - burnin);

    size_t samples = values.size() - burnin;
    size_t maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);

    double* gammaStat = new double[maxLag];
    // setting values to 0
    for (size_t i=0; i<maxLag; i++) {
        gammaStat[i] = 0;
    }
    double varStat = 0.0;

    for (size_t lag = 0; lag < maxLag; lag++) {
        for (size_t j = 0; j < samples - lag; j++) {
            double del1 = values.at(burnin + j) - mean;
            double del2 = values.at(burnin + j + lag) - mean;
            gammaStat[lag] += (del1 * del2);
        }

        gammaStat[lag] /= ((double) (samples - lag));

        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            // fancy stopping criterion :)
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            }
            // stop
            else
                maxLag = lag;
        }
    }

    // effective sample size
    double ess = samples / (varStat / gammaStat[0]);

    delete[] gammaStat;
    return ess;
}
