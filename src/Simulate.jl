"""
simulate.jl - pathogen evolution and transmission dynamic simulation tools
Justin Angevaare
May 2015
"""

function create_population(init_seq::Vector{Nucleotide}, init_var::Array)
"""
Create an infection database.
`init_seq` is assigned to the "external" infection source.
Each column of the `init_var` is assigned to an individual
"""
  # exposure times, exposure source, infection times, recovery times, covariate times, sequence times
  events = Array[Array[[NaN], [NaN],  [NaN],  [NaN],  [NaN], [0.]]]

  # covariate history, sequence history
  history = Array[Array[[[fill(NaN, length(init_var[:,1]))]],[init_seq]]]

  # push individuals to these arrays.
  for r = 1:size(init_var,2)
    push!(events, Array[Vector{Float64}[],    Vector{Int64}[],     Vector{Float64}[],     Vector{Float64}[],     [0.],   Vector{Float64}[]])
    push!(history, Array[[[init_var[:,r]]], Array{Nucleotide}[]])
  end

  # save as a population object type
  return Population(events, history)
end

function create_powerlaw(α::Float64, β::Float64, γ::Float64, η::Float64, dist=Euclidean())
  """
  This function creates a full parameterized power law function
  """
  @assert(α > 0, "invalid α specification")
  @assert(β > 0, "invalid β specification")
  @assert(γ > 0, "invalid γ specification")
  @assert(η > 0, "invalid η specification")

  return function(population::Population, source::Int, target::Int)
    """
    This simple `susceptibility_fun` returns the rate parameter for a `target` individuals from `source` individuals using the power law kernel with parameters α and β. Location must be specified with matching but arbitrary dimensions for each individual; specifically, each individual is represented by a column in an array. Distance by default is Euclidean, but any of the distance calculations in the Distance.jl package may be used.

    A zero distance is assigned a rate of γ, and a NaN distance (external source of infection) is assigned a rate of η.

    It's important to note that this function does not check the disease status of any individuals.

    This function also serves as a model for any user defined
    """
    # Ensure source is infectious that target hasn't been previously exposed (SIR model)
    if length(population.events[source][3]) > length(population.events[source][4]) && length(population.events[target][1]) == 0
    # Ensure source is infectious that target is susceptible (SIS* model)
    #if length(population.events[source][3]) > length(population.events[source][4]) && length(population.events[target][1]) == length(population.events[target][4])
      distance = evaluate(dist, population.history[source][1], population.history[target][1])
      if distance == 0.
        return γ
      elseif isnan(distance)
        return η
      else
        return α*distance^-β
      end
    else
      return 0.
    end
  end
end

function create_constantrate(τ::Float64)
  """
  Creates generic constant rate function
  """
  @assert(τ > 0, "invalid τ specification")
  return function(population::Population, individual::Int)
    """
    Provides a constant rate for latency_fun or recovery_fun
    """
    return τ
  end
end

function create_ratearray(population::Population, susceptibility_fun::Function, substitution_matrix::Array)
  """
  Generate an array which contains rates (for exponential distribution) for movement from between disease states, and mutation.
  `susceptibility_fun` is a function which generates a rate for each pair of target
  `substitution_matrix` is a 4x4 array containing single nucleotide polymorphism substitution rates
  """
  # Set up an array of zeros with rows for each potential source of exposure, for infection, recovery, and mutation at each base location, and columns for each individual...
  rate_array = RateArray(fill(0., (length(population.events)+1+1+length(population.history[1][2]), length(population.events))),
                         fill((), (length(population.events)+1+1+length(population.history[1][2]), length(population.events))))

  # Define events
  for r in 1:size(rate_array.events, 1)
    for c in 1:size(rate_array.events, 2)
      if r <= length(population.events)
        # Exposure event (from susceptible to exposed state)
        rate_array.events[r,c] = (1,c,r)
      elseif r == length(population.events)+1
        # Symptom event (from exposed to infectious state)
        rate_array.events[r,c] = (2,c)
      elseif r == length(population.events)+1+1
        # Recovery event (from infectious to recovered or susceptible* state)
        rate_array.events[r,c] = (3,c)
      else
        # Mutation event
        rate_array.events[r,c] = (4,c,r-(length(population.events)+1+1))
      end
    end
  end

  # External exposure rate
  for i = 2:size(rate_array.rates,2)
    rate_array.rates[1,i] = susceptibility_fun(pop, 1, i)
  end

  # Mutation of external pathogen
  rate_ref = sum(substitution_matrix,2)[:]
  nucleotide_ref = nucleotide("AGCU")
  for i = 1:length(population.history[1][2])
    rate_array.rates[length(population.events)+1+1+i,1] = rate_ref[findfirst(population.history[1][2][i] .== nucleotide_ref)]
  end

  # Return RateArray
  return rate_array
end

function onestep!(rate_array::RateArray, population::Population, time::Float64, susceptibility_fun::Function, latency_fun::Function, recovery_fun::Function, substitution_matrix::Array)
  """
  One event occurs, and appropriate updates are made to the RateArray and Population
  """
  rate_total = [0, cumsum(rate_array.rates[:])]

  if isinf(rate_total[end])
    increment = 0.
    event = rate_array.events[findfirst(isinf(rate_total))-1]

  else
    increment = rand(Exponential(1/rate_total[end]))
    event = rate_array.events[findfirst(total .> rand()*total[end])]
  end

  if event[1] == 1
    # S => E
    # Update rates - clear all individual specific rates
    rate_array.rates[:,event[2]] = 0.
    # Update rates - latency
    rate_array.rates[size(rate_array.rates,2)+1, event[2]] = latency_fun(population, event[2])
    # Update population - exposure time
    push!(population.events[event[2]][2], time+increment)
    # Update population - exposure source
    push!(population.events[event[2]][3], event[3])
    # Update population - sequence
    push!(population.history[event[2]][2], population.history[event[3]][2][end])
    # Update population - sequence time
    push!(population.events[event[2]][6], time+increment)
    # Update rates - mutation rates
    rate_ref = sum(substitution_matrix,2)[:]
    nucleotide_ref = nucleotide("AGCU")
    for i = 1:length(population.history[event[2]][2])
      rate_array.rates[size(rate_array.rates,2)+1+1+i, event[2]] = rate_ref[findfirst(population.history[event[2]][2][i] .== nucleotide_ref)]
    end

  elseif event[1] == 2
    # E => I
    # Update rates - clear latency
    rate_array.rates[event[3], event[2]] = 0.
    # Update rates - recovery
    rate_array.rates[size(rate_array.rates,2)+2, event[2]] = recovery_fun(population, event[2])
    # Update population - infection time
    push!(population.events[event[2]][4], time+increment)
    # Update rates - susceptibilities
    for i in 2:size(rate_array.rates,2)
      rate_array.rates[event[2],i] = susceptibility_fun(population, event[2], i)
    end

  elseif event[1] == 3
    # I => S*
    # Update rates - clear all individual specific rates (mutation and recovery)
    rate_array.rates[:,event[2]] = 0.
    # Update rates - susceptibilites of all other individuals
    rate_array.rates[event[2],:] = 0.
    # Update population - recovery time
    push!(population.events[event[2]][5], time+increment)
    # Update susceptibility* (under SIR framework, susceptibilities of 0 will be generated)
    for i = 1:size(rate_array.rates,2)
      rate_array[i, event[2]] = susceptibility_fun(population, i, event[2])
    end

  else
    # Mutation
    # Update population - sequence
    push!(population.history[event[2]][2], population.history[event[2]][2][end])
    nucleotide_ref = nucleotide("AGCU")
    nucleotide_mutation = substitution_matrix[:,findfirst(population.history[event[2]][2][end][event[3]] .== nucleotide_ref)]
    population.history[event[2]][2][end][event[3]] = nucleotide_ref[findfirst(rand(Multinomial(1, nucleotide_mutation/sum(nucleotide_mutation)),1))]
    # Update population - sequence time
    push!(population.events[event[2]][6], time+increment)
    # Update rates - mutation rates
    rate_ref = sum(substitution_matrix,2)[:]
    rate_array.rates[size(rate_array.rates,2)+1+1+event[3], event[2]] = rate_ref[findfirst(population.history[event[2]][2][end][event[3]] .== nucleotide_ref)]
  end
end

