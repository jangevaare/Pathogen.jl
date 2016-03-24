# Pathogen.jl
[![Build Status](https://travis-ci.org/jangevaare/Pathogen.jl.svg?branch=alpha)](https://travis-ci.org/jangevaare/Pathogen.jl)


## Introduction

Pathogen.jl is a package that provides utilities for the simulation and inference of pathogen phylodynamics, built in the [Julia language](http://julialang.org). Specifically, Pathogen.jl presents an extension to the individual level infectious disease transmission models (ILMs) of Deardon et al. (2010), to simultaneously model infectious disease transmission and evolution. Pathogen genomic sequences are used in conjunction with the covariate and disease state information of individuals to infer disease transmission pathways, external disease pressure, and infectivity parameters.

Pathogen.jl utilizes the [PhyloTrees.jl]("https://github.com/jangevaare/PhyloTrees.jl") package for pathogen sequence simulation and calculation of phylogenetic tree likelihoods.


## Package installation


    Pkg.update()
    Pkg.clone("https://github.com/jangevaare/Pathogen.jl/tree/alpha")



## Package usage

1. Create a population data frame containing covariate information for each individual in a population. This data frame should contain the location of each individual. Each row in this data frame is considered to represent a unique individual.

       using DataFrames
       population = Data.Frame(x = x_coordinates,
                               y = y_coordinates,
                               age = age)

2. Generate the event data frame using the Pathogen.jl package.

       using Pathogen
       events = create_events(population)

3. Create the following risk functions:

       function susceptibility(Θ::Vector{Float64}, population::DataFrame, i::Int64)
         # function of risk factors associated with the susceptibility of a susceptible
         # individual i
       end

       function transmissibility(Θ::Vector{Float64}, population::DataFrame, k::Int64)
         # function of risk factors associated with transmission of an infection from an
         # infectious individual k
       end

       function infectivity(Θ::Vector{Float64}, population::DataFrame, i::Int64, k::Int64)
         # function of risk factors involving both a susceptible individual i and an
         # infectious individual k
       end

       function sparks(Θ::Vector{Float64}, population::DataFrame, i::Int64)
         # function of risk factors associated with transmission of an infection from an
         # external disease source to a susceptible individual i
       end

       function latency(Θ::Vector{Float64}, population::DataFrame, j::Int64)
         # function of risk factors associated with the latency of an infection in an
         # exposed individual j
       end

       function removal(Θ::Vector{Float64}, population::DataFrame, k::Int64)
         # function of risk factors associated with the removal of an infectious
         # individual k
       end
