#include "MCMC.h"
#include "Model.h"
#include "DataWriter.h"
#include "Log.h"


MCMC::MCMC(Model& model, int numberOfGenerations, DataWriter& dataWriter) :
    _model(model), _numberOfGenerations(numberOfGenerations),
    _dataWriter(dataWriter)
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
        // TODO: MCMC should decide whether to accept or reject proposal
        _model.proposeNewState();
        _dataWriter.writeData(generation);
    }
}
