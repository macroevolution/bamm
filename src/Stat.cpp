#include "Stat.h"
#include <vector>
#include <cmath>


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
