# Pathogen

[![Build Status](https://travis-ci.org/jangevaa/Pathogen.jl.svg?branch=master)](https://travis-ci.org/jangevaa/Pathogen.jl)

Pathogen.jl is a package that provides utilities for the simulation and
inference of pathogen phylodynamics, built in the [Julia
language](http://julialang.org). Specifically, Pathogen.jl presents an extension
to the individual level infectious diesase transmission models (ILMs) of Deardon
et al. (2010), to simultaneously model infectious disease transmission and
evolution. Pathogen genomic sequences are used in conjunction with the covariate
and disease state information of individuals to infer disease transmission
pathways, external disease pressure, and infectivity parameters.
