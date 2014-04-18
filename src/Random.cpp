#include "Random.h"


Random::Random()
{
    warmUp();
}


Random::Random(unsigned long int seed) : _rng((long int)seed)
{
    warmUp();
}


void Random::warmUp()
{
    // Random number generations often need to be "warmed up"
    // (i.e., first few numbers thrown away)
    for (int i = 0; i < 1000; i++)
        _rng.uniformRv();
}


void Random::setSeed(unsigned long int seed)
{
    _rng.setSeed((long int)seed);
}


unsigned long int Random::getSeed() const
{
    return (unsigned long int)_rng.getSeed();
}


double Random::uniform()
{
    return _rng.uniformRv();
}


double Random::uniform(double a, double b)
{
    return _rng.uniformRv(a, b);
}


// Returns a random number in [a, b] under a uniform distribution
int Random::uniformInteger(int a, int b)
{
    return _rng.discreteUniformRv(a, b);
}


double Random::normal(double mean, double sd)
{
    return _rng.normalRv(mean, sd);
}


double Random::exponential(double rate)
{
    return _rng.exponentialRv(rate);
}


bool Random::trueWithProbability(double p)
{
    return uniform() < p;
}
