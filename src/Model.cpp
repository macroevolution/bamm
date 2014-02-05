#include "Model.h"
#include "MbRandom.h"
#include "Tree.h"
#include "Settings.h"
#include "Prior.h"


double Model::mhColdness = 1.0;


Model::Model(MbRandom* rng, Tree* tree, Settings* settings, Prior* prior) :
    _rng(rng), _tree(tree), _settings(settings), _prior(prior), _gen(0)
{
    // Reduce weird autocorrelation of values at start by calling RNG
    // a few times. TODO: Why is there a weird autocorrelation?
    for (int i = 0; i < 100; i++)
        _rng->uniformRv();

    // Event location scale is relative to the maximum root-to-tip length
    _scale = _settings->getUpdateEventLocationScale() *
        _tree->maxRootToTipLength();

    _updateEventRateScale = _settings->getUpdateEventRateScale();
    _localGlobalMoveRatio = _settings->getLocalGlobalMoveRatio();
    
    _poissonRatePrior = _settings->getPoissonRatePrior();

    // Initialize event rate to generate expected number of prior events
    _eventRate = 1 / _settings->getPoissonRatePrior();

    _acceptCount = 0;
    _rejectCount = 0;
    _acceptLast = -1;

    _lastDeletedEventMapTime = 0;
}


Model::~Model()
{
}
