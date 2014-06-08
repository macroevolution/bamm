Time-flip proposal: Time-constant and time-variable shifts
==========================================================

BAMM 2.0 introduces the concept of the "time-flip" proposal.
Previously, all shift events in BAMM were time-variable
(i.e., vary through time).
A time-variable shift event models speciation rate as

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

One problem in estimating the speciation rate through time
is that the statistical power provided by a tree is lower
at the start of the tree than at the end of the tree.
This is because there are fewer branches at the start of a tree
than at the end.
As a result, the initial speciation rate (:math:`\lambda_0`)
is more free to vary than the speciation rate at the end of the tree.
The value of :math:`k` is forced to compensate
for an inaccurate :math:`\lambda_0`.
This problem is compounded if the prior on :math:`\lambda_0`
does not match the true initial speciation rate.

If the true speciation rate of the phylogeny being analyzed
is constant through time, :math:`\lambda_0` should be estimated
as the speciation rate and :math:`k` should be estimated as 0.
However, :math:`\lambda_0` is often estimated inaccurately
due to the problem described above.
As a result, :math:`k` will be forced to be different than 0
in order for the speciation rate at the end to be accurate.
Our tests show that there is some bias in the estimation
of :math:`\lambda_0` and :math:`k` for time-constant trees.

Therefore, we introduced a "time-flip" proposal,
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

    \text{min}\left\{ 1, \cfrac{f(\lambda_0, k)\pi(\lambda_0, k)}
        {f(\lambda)\pi(\lambda)}
        \left| \cfrac{\partial (\lambda_0, k)}{\partial (\lambda, u)} \right|
    \right\}

where :math:`f(\lambda)` and :math:`\pi(\lambda)`
are the posterior and prior probabilities with the time-contant event,
:math:`f(\lambda_0, k)` and :math:`\pi(\lambda_0, k)`
are the posterior and prior probabilities with the time-variable event, and
:math:`\left| \cfrac{\partial (\lambda_0, k)}{\partial (\lambda, u)} \right|`
is the determinant of the Jacobian matrix for the transition from a
time-constant event to a time-variable event.
The determinant of this Jacobian matrix is

.. math::

    \left| \cfrac{\partial (\lambda_0, k)}{\partial (\lambda, u)} \right| =
    \renewcommand\arraystretch{2.3}
    \left| \begin{array}{cc}
        \cfrac{\partial \lambda_0}{\partial \lambda} &
        \cfrac{\partial \lambda_0}{\partial u} \\
        \cfrac{\partial k}{\partial \lambda} &
        \cfrac{\partial k}{\partial u}
    \end{array} \right| =
    \renewcommand\arraystretch{1.6}
    \left| \begin{array}{cc}
        \cfrac{\partial \lambda_0}{\partial \lambda} &
        \cfrac{\partial \lambda_0}{\partial u} \\
        0 & 1
    \end{array} \right| =
    \cfrac{\partial \lambda_0}{\partial \lambda}

From the equation above for :math:`\lambda_0` and equating
:math:`\lambda = \bar{\lambda}`, we obtain

.. math::

    \cfrac{\partial \lambda_0}{\partial \lambda} =
    \renewcommand\arraystretch{2.3}
    \left\{ \begin{array}{lr}
        \cfrac{kT}{e^{kT} - 1}        & \text{if } k < 0 \\
        \cfrac{kT}{2kT + e^{-kT} - 1} & \text{if } k > 0 \\
        1                             & \text{if } k = 0
    \end{array}
    \right.

The acceptance probability for a time-flip proposal
from a time-variable event to a time-constant event is

.. math::

    \text{min}\left\{ 1, \cfrac{f(\lambda)\pi(\lambda)}
        {f(\lambda_0, k)\pi(\lambda_0, k)}
        \left| \cfrac{\partial(\lambda_0, k)}{\partial(\lambda, u)}\right|^{-1}
    \right\}

which can be computed in the same way as described above.

Time-flip proposal options
--------------------------

The frequency in which a time-flip proposal occurs,
relative to other proposals, is given by ``updateRateLambdaTimeMode``
and ``updateRateBetaTimeMode`` for the diversification
and phenotypic evolution models, respectively.
When a new event is added to the tree, the probability that it is a
time-variable event is defined by ``lambdaIsTimeVariablePrior``.
The root event is assumed to be time-constant if ``lambdaShift0`` is 0;
otherwise, it is time-variable.
