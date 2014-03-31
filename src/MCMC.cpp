#include "MCMC.h"
#include "MbRandom.h"
#include "Model.h"
#include "DataWriter.h"
#include "Log.h"


MCMC::MCMC(MbRandom& rng, Model& model, int numberOfGenerations,
    DataWriter& dataWriter) : _rng(rng), _model(model),
    _numberOfGenerations(numberOfGenerations), _dataWriter(dataWriter)
{
}


MCMC::~MCMC()
{
}


void MCMC::run()
{
    log() << "\nRunning MCMC chain for "
          << _numberOfGenerations << " generations.\n";

    for (int generation = 0; generation < _numberOfGenerations; generation++) {
        _model.proposeNewState();

        double acceptanceRatio = _model.acceptanceRatio();
        if (_rng.uniformRv() < acceptanceRatio) {
            _model.acceptProposal();
        } else {
            _model.rejectProposal();
        }

        _dataWriter.writeData(generation);
    }
}
