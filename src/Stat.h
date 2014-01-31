#ifndef STAT_H
#define STAT_H

#include <vector>


class Stat
{
public:

    static double standard_deviation(const std::vector<double>& values);
    static double variance(const std::vector<double>& values);
};


#endif
