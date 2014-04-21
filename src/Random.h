#ifndef RANDOM_H
#define RANDOM_H


#include "MbRandom.h"


class Random
{
public:

    Random();
    Random(unsigned long int seed);

    void setSeed(unsigned long int seed);
    unsigned long int getSeed() const;

    double uniform();
    double uniform(double a, double b);

    int uniformInteger(int a, int b);

    double normal(double mean, double sd);
    double exponential(double rate);

    bool trueWithProbability(double p);

private:

    void warmUp();

    // Mutable allows const methods to call methods in MbRandom
    mutable MbRandom _rng;

    unsigned long int _seed;
};


#endif
