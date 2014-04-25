#include "LambdaTimeModeProposal.h"

#include "Model.h"
#include "Tree.h"
#include "SpExBranchEvent.h"

#include <cmath>

class Random;
class Settings;


LambdaTimeModeProposal::LambdaTimeModeProposal
    (Random& random, Settings& settings, Model& model)
    : _model(model), _tree(_model.getTreePtr()), _prior(random, &settings)
{
}


void LambdaTimeModeProposal::propose()
{
    _event = static_cast<SpExBranchEvent*>(_model.chooseEventAtRandom(true));

    _currentLambdaInit = _event->getLamInit();
    _currentLambdaShift = _event->getLamShift();
    _currentIsTimeVariable = _event->isTimeVariable();

    _currentLogLikelihood = _model.getCurrentLogLikelihood();
    _currentLogPrior = _model.computeLogPrior();

    if (_event->isTimeVariable()) {
        makeTimeConstant(_event);
        _lastTimeModeProposal = TimeConstant;
    } else {
        makeTimeVariable(_event);
        _lastTimeModeProposal = TimeVariable;
    }

    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();

    _proposedLogLikelihood = _model.computeLogLikelihood();
    _proposedLogPrior = _model.computeLogPrior();
}


void LambdaTimeModeProposal::makeTimeConstant(SpExBranchEvent* event)
{
    double lambdaInit = event->getLamInit();
    double lambdaShift = event->getLamShift();
    double timeLength = event->getTimeToTip();

    double meanLambda =
        computeMeanLambda(lambdaInit, lambdaShift, timeLength);

    event->setLamInit(meanLambda);
    event->setLamShift(0.0);
    event->setTimeVariable(false);
}


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

double LambdaTimeModeProposal::computeMeanLambda
    (double lam0, double k, double T)
{
    if (k > 0) {
        return 2.0 * lam0 + (lam0 / (k * T)) * (std::exp(-k * T) - 1.0);
    } else if (k < 0) {
        return (lam0 / (k * T)) * (std::exp(k * T) - 1.0);
    } else {
        return lam0;
    }
}


void LambdaTimeModeProposal::makeTimeVariable(SpExBranchEvent* event)
{
    double oldLambdaInit = event->getLamInit();
    double newLambdaShift = _prior.generateLambdaShiftFromPrior();
    double timeLength = event->getTimeToTip();

    double newLambdaInit =
        computeLambdaInit(oldLambdaInit, newLambdaShift, timeLength);

    event->setLamInit(newLambdaInit);
    event->setLamShift(newLambdaShift);
    event->setTimeVariable(true);
}


// The mean lambda (see "computeMeanLambda" above) is
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

double LambdaTimeModeProposal::computeLambdaInit
    (double lam, double k, double T)
{
    if (k > 0) {
        return lam / (2.0 + (std::exp(-k * T) - 1.0) / (k * T));
    } else if (k < 0) {
        return lam * k * T / (std::exp(k * T) - 1.0);
    } else {
        return lam;
    }
}


void LambdaTimeModeProposal::accept()
{
    _model.setCurrentLogLikelihood(_proposedLogLikelihood);
}


void LambdaTimeModeProposal::reject()
{
    _event->setLamInit(_currentLambdaInit);
    _event->setLamShift(_currentLambdaShift);
    _event->setTimeVariable(_currentIsTimeVariable);

    _tree->setNodeSpeciationParameters();
    _tree->setNodeExtinctionParameters();
}


double LambdaTimeModeProposal::acceptanceRatio()
{
    double logLikelihoodRatio = computeLogLikelihoodRatio();
    double logPriorRatio = computeLogPriorRatio();
    double logJacobian = computeLogJacobian();

    double t = _model.getTemperatureMH();
    double logRatio = t * (logLikelihoodRatio + logPriorRatio) + logJacobian;

    return std::min(1.0, std::exp(logRatio));
}


double LambdaTimeModeProposal::computeLogLikelihoodRatio()
{
    return _proposedLogLikelihood - _currentLogLikelihood;
}


double LambdaTimeModeProposal::computeLogPriorRatio()
{
    return _proposedLogPrior - _currentLogPrior;
}


double LambdaTimeModeProposal::computeLogJacobian()
{
    double k = _event->getLamShift();
    double T = _event->getTimeToTip();

    double jacobian = computeJacobian(k, T);
    double logJacobian = std::log(jacobian);

    if (_lastTimeModeProposal == TimeConstant) {
        return logJacobian;
    } else {
        return -logJacobian;
    }
}


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
// The initial lambda (see "calculateLambdaInit" above) is:
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

double LambdaTimeModeProposal::computeJacobian(double k, double T)
{
    if (k > 0) {
        return 1.0 / ((2.0 + std::exp(-k * T) - 1.0) / (k * T));
    } else if (k < 0) {
        return (k * T) / (std::exp(k * T) - 1.0);
    } else {
        return 1.0;
    }
}
