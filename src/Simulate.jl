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
  events = Array[Array[[NaN], [NaN],  [NaN],  [NaN],  [NaN], [0]],
                 Array[[],    [],     [],     [],     [0],   []]]

  # covariate history, sequence history
  history = Array[Array[[[fill(NaN, length(init_var[:,1]))]],[init_seq]],
                  Array[[[init_var[:,1]]],[]]]

  # push individuals to these arrays.
  for r = 2:size(init_var,2)
    push!(events, Array[[], [], [], [], [0], []])
    push!(history, Array[[[init_var[:,r]]],[]])
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
    distance = evaluate(dist, population.history[source][1], population.history[target][1])
    if distance == 0
      return γ
    elseif isnan(distance)
      return η
    else
      return α*distance^-β
    end
  end
end

function create_constantrate(τ::Float64)
  """
  Creates generic constant rate function
  """
  @assert(τ > 0, "invalid τ specification")
  return function()
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

  return rate_array
end

function onestep!(rates::RateArray, pop::Population, time::Float, susceptibility_fun::Function, latency_fun::Function, recovery_fun::Function, substitution_matrix::Array)
  """
  One event occurs, and appropriate updates are made to the RateArray and Population
  """
  increment = rand(Exponential(1/sum(rates.rates)))
  total = [0, cumsum(rates.rates[:])]
  event = rates.events[findfirst(total .> rand()*total[end])]
end
