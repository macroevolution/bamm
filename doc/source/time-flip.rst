.. _timefliptheory:

Time-flip proposal: Time-constant and time-variable shifts
==========================================================

**Note: The time-flip proposal is still being tested for reliability,
so it is currently turned off in the template and example control files.
It is available but is not the default. It is currently under development.**

BAMM 2.0 introduces the concept of the "time-flip" proposal.
Previously, the default model in BAMM assumed that all shift events in BAMM led to a distinct time-varying speciation (or phenotypic evolution) process. A time-variable shift event models speciation rate as

.. math::

    \renewcommand\arraystretch{1.3}
    \lambda(t) = \left\{
        \begin{array}{lr}
            \lambda_0 e^{kt}        & \text{if } k < 0   \\
            \lambda_0 (2 - e^{-kt}) & \text{if } k > 0 \\
            \lambda_0               & \text{if } k = 0
        \end{array}
    \right.

where :math:`\lambda_0` is the initial speciation rate,
:math:`k` is the rate parameter for the speciation rate,
and :math:`t` is the time since the shift event began.

We have found that this default assumption can lead to biased inference on speciation (or phenotypic evolutionary rates) in some areas of parameter space. In BAMM 2.0, we do not force rate shift events to be time-varying (as in BAMM 1.0) or constant through time (as in MEDUSA). Rather, we allow diversification submodels to "flip" between time-constant and time-variable rate modes. With this formulation, rate shifts will lead to a constant-rate diversification process unless the data contains sufficient evidence for temporally varying macroevolutionary rate dynamics. A given rate partition in the data can toggle between time-varying and constant-rate models in proportion to the posterior probability that the true process includes rate variation through time. We have found that this model substantially improves the performance of BAMM. 

The mechanics of the updated BAMM model entail the concept of a "time-flip" proposal,
which flips the time mode of a randomly-chosen shift event
from/to time-variable and time-constant modes.
A time-constant shift event models speciation as

.. math::

    \lambda(t) = \lambda

where :math:`\lambda` is the constant speciation rate.
With this model, the bias in the initial speciation rate is reduced
because it is not controlled by a separate parameter (:math:`\lambda_0`).
Instead, the entire speciation rate across the tree is controlled
by a single parameter (:math:`\lambda`).

Because the rate function is the same for both the diversification model
(:math:`\lambda`) and the phenotypic evolution model (:math:`\beta`),
this discussion will focus on the diversification model.


Calculation of parameters after a time-flip
-------------------------------------------

When flipping between time-constant and time-variable modes,
the rate function parameters are modified such that
the mean speciation rate through time is the same for both modes
(the time starts at the event location and ends at the tip of the tree).

If the chosen shift event is time-variable and the proposal to flip
the mode to time-constant is accepted, then the new speciation rate
is set to the mean value of the old speciation rate
(calculated from the time the event started until the tip of the tree):

.. math::

    \bar{\lambda} = \frac{1}{T} \int\limits_0^T \lambda (t) dt

where :math:`T` is the time elapsed since the beginning of the shift event.
The solution to this integral is

.. math::

    \renewcommand\arraystretch{2.3}
    \bar{\lambda} = \left\{
        \begin{array}{lr}
            \cfrac{\lambda _{0}}{kT} (e^{kT} - 1)        & \text{if } k < 0 \\
            \cfrac{\lambda _{0}}{kT} (2kT + e^{-kT} - 1) & \text{if } k > 0 \\
            \lambda _{0}                                 & \text{if } k = 0
        \end{array}
    \right.

If the chosen shift event is time-constant and the proposal to flip
the mode to time-variable is accepted, then the :math:`k` parameter is chosen
from its prior distribution. The :math:`\lambda _{0}` parameter
is calculated from the mean :math:`\lambda` (above equation).
Solving for :math:`\lambda _{0}` gives

.. math::

    \renewcommand\arraystretch{2.3}
    \lambda _{0} = \left\{
        \begin{array}{lr}
            \cfrac{\bar{\lambda}kT}{e^{kT} - 1}        & \text{if } k < 0 \\
            \cfrac{\bar{\lambda}kT}{2kT + e^{-kT} - 1} & \text{if } k > 0 \\
            \bar{\lambda}                              & \text{if } k = 0
        \end{array}
    \right.

where :math:`\bar{\lambda}` is taken from the speciation rate
for the time-constant event.


Acceptance probability for a time-flip proposal
-----------------------------------------------

The acceptance probability for a time-flip proposal
from a time-constant event to a time-variable event is

.. math::

    \text{min}\left\{ 1, \cfrac{f(\lambda_0, k)}{f(\lambda)} \times
        \cfrac{\pi(\lambda_0, k)}{\pi(\lambda)} \times
        \cfrac{q(\lambda | \lambda_0, k)}{q(\lambda_0, k | \lambda)q(u)} \times
        \left| \cfrac{\partial g(\lambda, u)}{\partial (\lambda, u)} \right|
    \right\}

where :math:`f(\lambda)` and :math:`\pi(\lambda)`
are the posterior and prior probabilities for the time-contant event,
:math:`f(\lambda_0, k)` and :math:`\pi(\lambda_0, k)`
are the posterior and prior probabilities for the time-variable event,
:math:`q(\lambda | \lambda_0, k)` is the probability of proposing
a move to parameter :math:`\lambda` given that the current
parameters are :math:`\lambda_0` and :math:`k`,
:math:`q(\lambda_0, k | \lambda)` is the probability of the reverse move,
:math:`q(u)` is the probability density of the random variable :math:`u`, and
:math:`\left| \cfrac{\partial g(\lambda, u)}{\partial (\lambda, u)} \right|`
is the determinant of the Jacobian matrix for the transition from a
time-constant event to a time-variable event,
using the mapping function :math:`g`:

.. math::

    \left[ \begin{array}{c}
        \lambda_0 \\
        k
    \end{array} \right] =
    g(\lambda, u) =
    \left[ \begin{array}{c}
        g_1(\lambda, u) \\
        g_2(\lambda, u)
    \end{array} \right]

where

.. math::

    g_1(\lambda, u) =
    \renewcommand\arraystretch{2.3}
    \left\{ \begin{array}{lr}
        \cfrac{\lambda uT}{e^{uT} - 1}        & \text{if } u < 0 \\
        \cfrac{\lambda uT}{2uT + e^{-uT} - 1} & \text{if } u > 0 \\
        \lambda                               & \text{if } u = 0
    \end{array} \right.

and :math:`g_2(\lambda, u) = u`.
The value :math:`u` is a random number taken from the prior distribution
of :math:`k`.
The determinant of the Jacobian matrix is therefore

.. math::

    \left| \cfrac{\partial g(\lambda, u)}{\partial (\lambda, u)} \right| =
    \renewcommand\arraystretch{2.3}
    \left| \begin{array}{cc}
        \cfrac{\partial g_1(\lambda, u)}{\partial \lambda} &
        \cfrac{\partial g_1(\lambda, u)}{\partial u} \\
        \cfrac{\partial g_2(\lambda, u)}{\partial \lambda} &
        \cfrac{\partial g_2(\lambda, u)}{\partial u}
    \end{array} \right| =
    \renewcommand\arraystretch{1.6}
    \left| \begin{array}{cc}
        \cfrac{\partial g_1(\lambda, u)}{\partial \lambda} &
        \cfrac{\partial g_1(\lambda, u)}{\partial u} \\
        0 & 1
    \end{array} \right| =
    \cfrac{\partial g_1(\lambda, u)}{\partial \lambda}

This partial derivative is easy to calculate:

.. math::

    \cfrac{\partial g_1(\lambda, u)}{\partial \lambda} =
    \renewcommand\arraystretch{2.3}
    \left\{ \begin{array}{lr}
        \cfrac{uT}{e^{uT} - 1}        & \text{if } u < 0 \\
        \cfrac{uT}{2uT + e^{-uT} - 1} & \text{if } u > 0 \\
        1                             & \text{if } u = 0
    \end{array}
    \right.

The acceptance probability for a time-flip proposal
from a time-variable event to a time-constant event is

.. math::

    \text{min}\left\{ 1,
        \cfrac{f(\lambda)}{f(\lambda_0, k)} \times
        \cfrac{\pi(\lambda)}{\pi(\lambda_0, k)} \times
        \cfrac{q(\lambda_0, k | \lambda)q(u)} {q(\lambda | \lambda_0, k)} \times
        \left|\cfrac{\partial g(\lambda_0,u)}{\partial(\lambda,u)}\right|^{-1}
    \right\}

Time-flip proposal options
--------------------------

The frequency in which a time-flip proposal occurs,
relative to other proposals, is given by ``updateRateLambdaTimeMode``
and ``updateRateBetaTimeMode`` for the diversification
and phenotypic evolution models, respectively.
When a new event is added to the tree, the probability that it is time-variable
is defined by ``lambdaIsTimeVariablePrior`` (or ``betaIsTimeVariablePrior``).
If the probability that an event is time-variable is between 0 and 1,
the initial root event is assumed to be time-constant if ``lambdaShift0`` is 0;
otherwise, it is time-variable.
A similar assumption is made for ``betaShiftInit``.

To constrain BAMM such that all diversification shifts lead to
time-varying processes only, set::

    lambdaIsTimeVariablePrior = 1
    updateRateLambdaTimeMode = 0

To constrain BAMM such that all diversification shifts lead to time-constant
diversification processes only, set::

    lambdaIsTimeVariablePrior = 0
    updateRateLambdaTimeMode = 0
    lambdaShift0 = 0

Make similar adjustments to the corresponding *beta* options for
the phenotypic evolution model type.
