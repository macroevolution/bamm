:orphan:

.. _mc3:

.. |MC3| replace:: (MC)\ :sup:`3`

Metropolis coupled Markov chain Monte Carlo [|MC3|]
===================================================

The Markov chain Monte Carlo (MCMC) method implemented in BAMM
estimates the posterior probability distribution of diversification models that best describe a particular phylogenetic dataset. Before version 2.0.0, BAMM used a single Markov chain
to explore the landscape of models and their parameter values.
A single chain, however, may get stuck in local optima,
which results in less mixing and more time needed for convergence.

Metropolis coupled Markov chain Monte Carlo [|MC3|],
introduced in BAMM 2.0.0, increases the number of Markov chains
that explore the landscape of models and their parameters.
Not only are there more chains, but each one has a different *temperature*.
The *cold* chain is the main chain and behaves as before.
The *heated* chains explore a landscape that is flatter than
the landscape explored by the cold chain.
Therefore, it is easier for a heated chain to cross deep valleys
in the landscape and not get stuck in local optima.

After the chains independently explore the landscape
for a certain number of generations,
a proposal is made to swap the temperatures of two random chains.
If a swap involving the cold chain is successful,
it will cause a previously heated chain to become the main chain.
As a result, a chain that is stuck in a local optimum
may immediately jump to another area of the landscape.

The implementation of |MC3| in BAMM follows that described in
`Altekar et al. 2004
<http://bioinformatics.oxfordjournals.org/content/20/3/407.full.pdf>`_.
For each chain :math:`i`, its temperature :math:`\beta_i` is set to
:math:`\beta_i = [1 + \Delta T \times (i - 1)]^{-1}`,
where :math:`\Delta T` is the temperature increment parameter.
For example, if there are 4 chains and :math:`\Delta T = 0.1`,
the temperatures of each chain are 1, 0.9091, 0.8333, and 0.7692.
The temperature of the cold chain is always 1.
Note that the value of :math:`\Delta T` should be greater than 0
and chosen such that the probability of accepting a swap
is between 20% and 60% (Altekar et al. 2004).

The temperature of each chain :math:`i` goes into the calculation
of the acceptance probability for a within-model proposal
(i.e., not involving changes in the dimensionality of the model):

.. math::

    \alpha_i = \text{min}\left\{ 1,
        \left(
        \cfrac{f(\theta_i')}{f(\theta_i)} \times
        \cfrac{\pi(\theta_i')}{\pi(\theta_i)}
        \right)^{\beta_i} \times
        \cfrac{q(\theta_i | \theta_i')}{q(\theta_i' | \theta_i)}
    \right\}

where :math:`\theta_i` and :math:`\theta_i'` are parameter vectors
corresponding to the current and proposed states for chain :math:`i`,
:math:`f` and :math:`\pi` are the corresponding likelihood
and prior density functions,
and :math:`q(\theta_i' | \theta_i)` is the relative probability
of proposing a move to parameter vector :math:`\theta_i'`
given that the current state is :math:`\theta_i`.
A similar calculation is done for the acceptance probability for proposals
that change the dimensionality of the model.

After a certain number of generations, two randomly chosen chains
:math:`j` and :math:`k` are swapped with acceptance probability

.. math::

    \alpha = \text{min}\left\{ 1,
        \left(\cfrac{f(\theta_k)}{f(\theta_j)}\right)^{\beta_j} \times
        \left(\cfrac{f(\theta_j)}{f(\theta_k)}\right)^{\beta_k}
    \right\}

|MC3| settings in BAMM
----------------------

The number of chains to use is specified with ``numberOfChains``.
The :math:`\Delta T` value is specified with ``deltaT``.
The number of generations between chain swap proposals
is specified with ``swapPeriod``.
BAMM will output to a file named in ``chainSwapFileName``
(the default is ``chain_swap.txt``)
the generation in which a swap proposal occurred,
the ranking of the two chains chosen
(1 is the cold chain, 2 is the first heated chain, etc.),
and whether the swap was accepted.
