# Pathogen.jl
[![Build Status](https://travis-ci.org/jangevaare/Pathogen.jl.svg?branch=master)](https://travis-ci.org/jangevaare/Pathogen.jl)

## Introduction

Pathogen.jl is a package that provides utilities for the simulation and inference of pathogen phylodynamics, built in the [Julia language](http://julialang.org). Specifically, Pathogen.jl presents an extension to the individual level infectious disease transmission models (ILMs) of Deardon et al. (2010), to simultaneously model infectious disease transmission and evolution. Pathogen genomic sequences are used in conjunction with the covariate and disease state information of individuals to infer disease transmission pathways, external disease pressure, and infectivity parameters.

Pathogen.jl utilizes the packages [PhyloTrees.jl](https://github.com/jangevaare/PhyloTrees.jl) and [PhyloModels.jl](https://github.com/jangevaare/PhyloModels.jl) for pathogen sequence simulation and calculation of phylogenetic tree likelihoods.


![Phylodynamic simulation](epianimation.gif?raw=true)


## Package installation


    Pkg.update()
    Pkg.clone("https://github.com/jangevaare/Pathogen.jl")
