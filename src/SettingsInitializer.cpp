#include "SettingsInitializer.h"

#include "Random.h"
#include "Settings.h"
#include "Node.h"
#include "Prior.h"
#include "SpExModel.h"
#include "TraitModel.h"
#include "Log.h"

#include <cstdlib>


SettingsInitializer::SettingsInitializer(Settings& settings) :
    _random(), _settings(settings), _tree(_random, _settings)
{
}


void SettingsInitializer::initializeSpExSettings()
{
    initSpExPriors();    // this must be called before initSpExInitParams()
    initSpExInitParams();
}


void SettingsInitializer::initSpExPriors()
{
    double TMax = _tree.maxRootToTipLength();
    double TMaxInverse = 1.0 / TMax;

    double nTaxa = _tree.getNumberTips();
    double sppProb = TMaxInverse * std::log(nTaxa / 2.0);

    double lambdaInitPrior = 1.0 / (sppProb * 5.0);
    double lambdaShiftPrior = -TMaxInverse * std::log(0.1) / 2.0;

    // Note: muInitPrior is correctly set to lambdaInitPrior
    _settings.setIfNotUserDefined("lambdaInitPrior", lambdaInitPrior);
    _settings.setIfNotUserDefined("lambdaShiftPrior", lambdaShiftPrior);
    _settings.setIfNotUserDefined("muInitPrior", lambdaInitPrior);
    _settings.setIfNotUserDefined("poissonRatePrior", 1.0);
}


void SettingsInitializer::initSpExInitParams()
{
    Prior prior(_random, &_settings);

    double logLikelihood = 0.0;
    int count = 0;

    do {
        double lambdaInit0 = prior.generateLambdaInitFromPrior();
        double muInit0 = _random.uniform() * lambdaInit0;

        _settings.setIfNotUserDefined("lambdaInit0", lambdaInit0);
        _settings.setIfNotUserDefined("muInit0", muInit0);

        SpExModel model(_random, _settings);
        logLikelihood = model.getCurrentLogLikelihood();
    } while (std::isinf(logLikelihood) && count++ < 1000);

    if (std::isinf(logLikelihood)) {
        log(Error) << "Could not find good initial parameters.\n";
        std::exit(1);
    }
}


void SettingsInitializer::initializeTraitSettings()
{
    initTraitPriors();    // this must be called before initSpExInitParams()
    initTraitParams();
}


void SettingsInitializer::initTraitPriors()
{
    double TMax = _tree.maxRootToTipLength();
    double TMaxInverse = 1.0 / TMax;

    double sigma2 = calculateSigma2();

    double betaInitPrior = 1.0 / (sigma2 * 5.0);
    double betaShiftPrior = -TMaxInverse * std::log(0.1) / 2.0;

    _settings.setIfNotUserDefined("betaInitPrior", betaInitPrior);
    _settings.setIfNotUserDefined("betaShiftPrior", betaShiftPrior);
}


// From O'Meara et al 2006 (note, sigma^2 is beta):
//
//               (X - a)' * C^-1 * (X - a)
//     sigma^2 = -------------------------
//                          N
//
// where
//
//     N is the number of taxa
//     X is a vector of traits (size is N)
//     C is an N x N matrix where (i, j) element is the time
//        from root to the MCRA of taxa i and j
//     1 is a vector of ones (size is N)
//     a = (1' * C^-1 * 1)^-1 * (1' * C^-1 * X)
//     Note: M' is the transpose of matrix M
//     Note: M^-1 is the inverse of matrix M

double SettingsInitializer::calculateSigma2()
{
    try {    // matrix operations could fail (e.g., inverse)
        double N = _tree.getNumberTips();
        vec X(_tree.traitValues());
        const mat& C = createC();
        mat CInv = C.i();    // inverse
        vec one(N, fill::ones);
        rowvec oneT = one.t();    // transpose

        mat aMat = (oneT * CInv * one).i() * (oneT * CInv * X);
        double a = aMat(0, 0);    // aMat should have only one element

        mat sigma2Mat = ((X - a).t() * CInv * (X - a)) / N;
        double sigma2 = sigma2Mat(0, 0);    // sigma2Mat should have one element

        return sigma2;
    } catch (...) {
        log(Error) << "Auto-initialization of priors failed. "
            "Set priors manually.\n";
        std::exit(1);
    }
}


// C is an N x N matrix where (i, j) element is the time
// from root to the MCRA of taxa i and j

mat SettingsInitializer::createC()
{
    const std::vector<Node*>& tips = _tree.terminalNodes();
    int N = tips.size();

    mat C(N, N);

    for (int i = 0; i < N; i++) {
        Node* tip_i = tips[i];
        for (int j = 0; j < N; j++) {
            Node* tip_j = tips[j];
            Node* MRCA = _tree.getNodeMRCA(tip_i->getName(), tip_j->getName());
            C(i, j) = MRCA->pathLengthToRoot();
        }
    }

    return C;
}


void SettingsInitializer::initTraitParams()
{
    Prior prior(_random, &_settings);

    double logLikelihood = 0.0;
    int count = 0;

    do {
        double betaInit = prior.generateBetaInitFromPrior();

        _settings.setIfNotUserDefined("betaInit", betaInit);

        TraitModel model(_random, _settings);
        logLikelihood = model.getCurrentLogLikelihood();
    } while (std::isinf(logLikelihood) && count++ < 1000);

    if (std::isinf(logLikelihood)) {
        log(Error) << "Could not find good initial parameters.\n";
        std::exit(1);
    }
}
