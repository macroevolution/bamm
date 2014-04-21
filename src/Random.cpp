#include "Random.h"


// For the default constructor, initialize the MbRandom object first
// for it to generate a seed based on the clock, then assign it internally
Random::Random() : _rng(), _seed(_rng.getSeed())
{
    warmUp();
}


Random::Random(unsigned long int seed) : _rng((long int)seed), _seed(seed)
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


unsigned long int Random::getSeed() const
{
    return _seed;
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
