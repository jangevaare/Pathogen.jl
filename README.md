# Pathogen.jl
[![DOI](https://zenodo.org/badge/35234698.svg)](https://zenodo.org/badge/latestdoi/35234698)
[![Latest Release](https://img.shields.io/github/release/jangevaare/Pathogen.jl.svg)](https://github.com/jangevaare/Pathogen.jl/releases/latest)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jangevaare/Pathogen.jl/blob/master/LICENSE)

[![test-lts](https://github.com/jangevaare/Pathogen.jl/actions/workflows/test-lts.yml/badge.svg)](https://github.com/jangevaare/Pathogen.jl/actions/workflows/test-lts.yml)
[![test-stable](https://github.com/jangevaare/Pathogen.jl/actions/workflows/test-stable.yml/badge.svg)](https://github.com/jangevaare/Pathogen.jl/actions/workflows/test-stable.yml)
[![test-nightly](https://github.com/jangevaare/Pathogen.jl/actions/workflows/test-nightly.yml/badge.svg)](https://github.com/jangevaare/Pathogen.jl/actions/workflows/test-nightly.yml)
[![codecov.io](http://codecov.io/github/jangevaare/Pathogen.jl/coverage.svg?branch=master)](http://codecov.io/github/jangevaare/Pathogen.jl?branch=master)

Authors: Justin Angevaare, Zeny Feng, Rob Deardon

![Epidemic curve](https://github.com/jangevaare/Pathogen.jl/raw/master/examples/SIR%20Simulation/epiplot.png)

Pathogen.jl is a Julia software package for individual level models of infectious diseases (Deardon et al, 2010). It's capabilities include stochastic simulation and Bayesian inference of SEIR, SEI, SIR, and SI individual level models, with fully customizable functions describing individual specific transition rates between disease states (i.e. form of, and relevant risk factors to, susceptibility, transmissibility, latency, removal, and sparks functions). Pathogen.jl is written purely in Julia, which enables this generality without incurring performance costs.

![MCMC](https://github.com/jangevaare/Pathogen.jl/raw/master/examples/SIR%20Simulation/posterior.png)

Pathogen.jl infers transmission pathways (i.e. who-infected-who). This inference is completed using a Gibbs step in our specialized MCMC algorithm. This specialized MCMC algorithm also performs event time data augmentation. A detailed overview of this algorithm can be found [here](https://arxiv.org/abs/2002.05850).

## Installation
The current release can be installed from the Julia REPL with:

```julia
pkg> add Pathogen
```

The development version (master branch) can be installed with:

```julia
pkg> add Pathogen#master
```

![Posterior Transmission Network](https://github.com/jangevaare/Pathogen.jl/raw/master/examples/SIR%20Simulation/posterior_tn.png)

Examples of Pathogen.jl workflow are included in the examples directory as a Jupyter notebooks.
1. [SIR simulation, inference, and visualization](https://github.com/jangevaare/Pathogen.jl/blob/master/examples/SIR%20Simulation/SIR%20TN-ILM%20Simulation%20and%20Inference.ipynb)
2. [Analysis of a Tomato Spotted Wilt Virus experimental epidemic](https://github.com/jangevaare/Pathogen.jl/blob/master/examples/Tomato%20Spotted%20Wilt%20Virus/TSWV.ipynb)
3. [Analysis of 1861 Hagelloch Measles outbreak](https://github.com/jangevaare/Pathogen.jl/blob/master/examples/1861%20Hagelloch%20Measles/1861%20Hagelloch.ipynb)
![Epidemic simulation](https://github.com/jangevaare/Pathogen.jl/raw/master/examples/SIR%20Simulation/epianimation.gif)

## Citation and more information
This package is detailed in [Pathogen.jl: Infectious Disease Transmission Network Modeling with Julia](https://doi.org/10.18637/jss.v104.i04), in the Journal of Statistical Software.

    @article{pathogenjl,
      title   = {Pathogen.jl: Infectious Disease Transmission Network Modeling with Julia},
      author  = {Angevaare, Justin and
                 Feng, Zeny and
                 Deardon, Rob},
      year    = {2022},
      journal = {Journal of Statistical Software},
      volume  = {104},
      number  = {4},
      pages   = {1–30},
      url     = {https://www.jstatsoft.org/index.php/jss/article/view/v104i04},
      doi     = {10.18637/jss.v104.i04}}

## Random number generators

Random simulation and inference methods accept an optional first positional
`rng::AbstractRNG`. Calls without it use Julia's default RNG. For example:

```julia
using Random
rng = MersenneTwister(123)
simulate!(rng, sim; tmax=100.0)
obs = observe(rng, sim, Uniform(); seq_len=100) # SI/SEI example
start!(rng, mcmc; attempts=100)
iterate!(rng, mcmc, 100, 1.0)
```

The supplied RNG controls event simulation, observation delays, pathogen
sequences, chain initialization, proposals, and acceptance decisions. The same
RNG-first convention applies to stochastic `generate` methods and the internal
`initialize`, `update!`, and `_accept` methods.

For multiple chains, `start!` and `iterate!` draw one `UInt64` seed per chain
from the supplied RNG in chain order and create independent `MersenneTwister`
streams on the workers. Results therefore do not depend on worker scheduling.
A single chain uses the supplied RNG directly. Keep and reuse the RNG when
continuing a run; chains do not store RNG state. Exact streams may differ across
Julia and dependency versions, and between single-chain and multi-chain runs.

## Inference performance

Inference reuses a likelihood workspace and two transmission-rate cache buffers
within each chain. Removal events subtract the departing infector's contribution
from susceptible hazards, and parameter proposals reuse event ordering when the
event times are unchanged. Dense pairwise risk functions still require quadratic
work and rate storage; these changes do not assume a distance cutoff or discard
small transmission risks.

Sampled `TransmissionNetwork.internal` matrices now use Julia's sparse Boolean
matrix storage. Matrix indexing and construction from dense Boolean matrices
remain supported. Posterior network distributions remain dense because different
samples can assign many possible parents to each individual.

The internal `loglikelihood` function accepts a `workspace` keyword. Outputs that
use its rate and transmission buffers are borrowed until the next evaluation
with that workspace; calls without a workspace return independent buffers.
Use a separate workspace for each concurrent evaluation. The legacy
`early_decision_value` keyword is accepted but does not prune partial likelihoods:
continuous-time log-density contributions can be positive.

Run `benchmark/likelihood.jl` in an environment containing this checkout to measure
warmed likelihood time and allocations at increasing population sizes. These
measurements exclude compilation and fixture construction.

