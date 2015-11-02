#include "TimeModeProposal.h"

#include "Model.h"
#include "Tree.h"
#include "BranchEvent.h"

#include <algorithm>
#include <cmath>

class Random;
class Settings;


TimeModeProposal::TimeModeProposal
    (Random& random, Settings& settings, Model& model)
    : _model(model), _tree(_model.getTreePtr()), _prior(random, &settings)
{
}


void TimeModeProposal::propose()
{
    _event = _model.chooseEventAtRandom(true);

    _currentInitParam = initialParameter(_event);
    _currentRateParam = rateParameter(_event);
    _currentIsTimeVariable = isTimeVariable(_event);

    _currentLogLikelihood = _model.getCurrentLogLikelihood();
    _currentLogPrior = _model.computeLogPrior();

    if (_currentIsTimeVariable) {
        makeTimeConstant(_event);
        _lastTimeModeProposal = TimeConstant;
    } else {
        makeTimeVariable(_event);
        _lastTimeModeProposal = TimeVariable;
    }

    setModelParameters();

    _proposedLogLikelihood = _model.computeLogLikelihood();
    _proposedLogPrior = _model.computeLogPrior();
}


void TimeModeProposal::makeTimeConstant(BranchEvent* event)
{
    double initParam = initialParameter(event);
    double rateParam = rateParameter(event);
    double timeLength = event->getTimeToTip();

    double meanRate = computeMeanRate(initParam, rateParam, timeLength);

    setEventParameters(event, meanRate, 0.0, false);
}


// The following is specific to lambda, but applies similarly to beta.
// The equation for lambda through time, lam(t), is
//
//                 /  lam0 (2 - e^(-kt))    k > 0
//                |
//     lam(t) =  <    lam0 e^(kt)           k < 0
//                |
//                 \  lam0                  k = 0
//
// where lam0 is the initial lambda, k is the rate parameter, and t is time.
// The mean lambda, named "lam" here, is calculated as
//
//              T
//           1 /
//     lam = - | lam(t) dt
//           T /
//              0
//
// where T is the time since the event to the tip of the tree.
// Solving this integral results in
//
//             /           lam0
//            |   2 lam0 + ---- (e^(-kT) - 1)    k > 0
//            |             kT
//            |
//     lam = <    lam0
//            |   ---- (e^(kT) - 1)              k < 0
//            |    kT
//            |
//             \  lam0                           k = 0

double TimeModeProposal::computeMeanRate(double init, double k, double T)
{
    if (k > 0) {
        return 2.0 * init + (init / (k * T)) * (std::exp(-k * T) - 1.0);
    } else if (k < 0) {
        return (init / (k * T)) * (std::exp(k * T) - 1.0);
    } else {
        return init;
    }
}


void TimeModeProposal::makeTimeVariable(BranchEvent* event)
{
    double oldInitParam = initialParameter(event);
    double newRateParam = rateParameterFromPrior();
    double timeLength = event->getTimeToTip();

    double newInitParam =
        computeRateInit(oldInitParam, newRateParam, timeLength);

    setEventParameters(event, newInitParam, newRateParam, true);
}


// The following is specific to lambda, but applies similarly to beta.
// The mean lambda (see "computeMeanRate" above) is
//
//             /           lam0
//            |   2 lam0 + ---- (e^(-kT) - 1)    k > 0
//            |             kT
//            |
//     lam = <    lam0
//            |   ---- (e^(kT) - 1)              k < 0
//            |    kT
//            |
//             \  lam0                           k = 0
//
// where lam0 is the initial lambda, k is the rate parameter,
// and T is time since the event to the tip of the tree.
//
// Solving for lam0 results in
//
//               /        lam
//              |   ---------------
//              |   2 + e^(-kT) - 1    k > 0
//              |       -----------
//              |           kT
//      lam0 = <
//              |     lam kT
//              |   ----------         k < 0
//              |   e^(kT) - 1
//              |
//               \  lam                k = 0

double TimeModeProposal::computeRateInit(double mean, double k, double T)
{
    if (k > 0) {
        return mean / (2.0 + (std::exp(-k * T) - 1.0) / (k * T));
    } else if (k < 0) {
        return mean * k * T / (std::exp(k * T) - 1.0);
    } else {
        return mean;
    }
}


void TimeModeProposal::accept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void TimeModeProposal::reject()
{
    setEventParameters(_event, _currentInitParam, _currentRateParam,
        _currentIsTimeVariable);

    setModelParameters();
}


double TimeModeProposal::acceptanceRatio()
{
    double logLikelihoodRatio = computeLogLikelihoodRatio();
    double logPriorRatio = computeLogPriorRatio();
    double logJacobian = computeLogJacobian();

    double t = _model.getTemperatureMH();
    double logRatio = t * (logLikelihoodRatio + logPriorRatio) + logJacobian;

    if (std::isfinite(logRatio)) {
        return std::min(1.0, std::exp(logRatio));
    } else {
        return 0.0;
    }
}


double TimeModeProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double TimeModeProposal::computeLogPriorRatio()
{
    return _proposedLogPrior - _currentLogPrior;
}


double TimeModeProposal::computeLogJacobian()
{
    double k = rateParameter(_event);
    double T = _event->getTimeToTip();

    double jacobian = computeJacobian(k, T);
    double logJacobian = std::log(jacobian);

    if (_lastTimeModeProposal == TimeConstant) {
        return logJacobian;
    } else {
        return -logJacobian;
    }
}


// The following is specific to lambda, but applies similarly to beta.
// The determinant of the Jacobian matrix for transitioning from a
// time constant to a time variable mode is
//
//           | d(lam0)  d(lam0) |
//           | -------  ------- |   | d(lam0)  d(lam0) |
//           |  d(lam)    d(u)  |   | -------  ------- |   d(lam0)
//     |J| = |                  | = |  d(lam)    d(u)  | = -------
//           |   d(k)     d(k)  |   |                  |    d(lam)
//           |  -----     ----  |   |    0         1   |
//           |  d(lam)    d(u)  |
//
// where d() is the partial derivative, lam0 is the initial lambda,
// k is the lambda shift, lam is the mean (constant) lambda, and
// u is a random variable from the same prior distribution as k.
//
// The initial lambda (see "computeRateInit" above) is:
//
//               /        lam
//              |   ---------------
//              |   2 + e^(-kT) - 1    k > 0
//              |       -----------
//              |           kT
//      lam0 = <
//              |     lam kT
//              |   ----------         k < 0
//              |   e^(kT) - 1
//              |
//               \  lam                k = 0
//
// Therefore,
//
//                  /         1
//                 |   ---------------
//                 |   2 + e^(-kT) - 1    k > 0
//                 |       -----------
//      d(lam0)    |           kT
//      ------- = <
//       d(lam)    |       kT
//                 |   ----------         k < 0
//                 |   e^(kT) - 1
//                 |
//                  \  1                  k = 0
//

double TimeModeProposal::computeJacobian(double k, double T)
{
    if (k > 0) {
        return 1.0 / ((2.0 + std::exp(-k * T) - 1.0) / (k * T));
    } else if (k < 0) {
        return (k * T) / (std::exp(k * T) - 1.0);
    } else {
        return 1.0;
    }
}
