#ifndef SETTINGS_INITIALIZER_H
#define SETTINGS_INITIALIZER_H


#include "Random.h"
#include "Tree.h"

#include <armadillo>
using namespace arma;

class Settings;


class SettingsInitializer
{
public:

    SettingsInitializer(Settings& settings);

    void initializeSpExSettings();
    void initializeTraitSettings();

private:

    void initSpExPriors();
    void initSpExInitParams();

    void initTraitPriors();
    double calculateSigma2();
    mat createC();
    void initTraitParams();

    Random _random;
    Settings& _settings;
    Tree _tree;
};


#endif
