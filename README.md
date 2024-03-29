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
