"""
Generate a nucleotide sequence of length `n`, with specific nucleotide frequencies
"""
function create_seq(n::Int, π_A::Float64, π_T::Float64, π_C::Float64, π_G::Float64)
  @assert(sum([π_A, π_T, π_C, π_G]) == 1, "Nucleotide frequencies must sum to 1")
  @assert(all(0 .< [π_A, π_T, π_C, π_G] .< 1), "Each nucleotide frequency must be between 0 and 1")
  return convert(Nucleotide2bitSeq, findn(rand(Multinomial(1, [π_A, π_T, π_C, π_G]),n))[1])
end


"""
Create an infection database.
`init_seq` is assigned to the "external" infection source.
Each column of the `init_var` is assigned to an individual
"""
function create_population(init_seq::Nucleotide2bitSeq, init_var::Array)
  # exposure times, exposure source, infection times, recovery times, covariate times, sequence times
  events = Array[Array[[NaN], [NaN],  [NaN],  [NaN],  [NaN], [0.]]]

  # covariate history, sequence history
  history = Array[Array[Array[fill(NaN, length(init_var[:,1]))],[init_seq]]]

  # event time, event type
  timeline = Array[[0.],[(0,0,0)]]

  # push individuals to these arrays.
  for r = 1:size(init_var,2)
    push!(events, Array[Float64[],    Int64[],     Float64[],     Float64[],     [0.],   Float64[]])
    push!(history, Array[Array[init_var[:,r]], Vector{Nucleotide2bitSeq}[]])
  end

  # save as a population object type
  return Population(events, history, timeline)
end


"""
This function creates a full parameterized power law function
"""
function create_powerlaw(α::Float64, β::Float64, η::Float64, dist=Euclidean()::Metric)
  @assert(α > 0, "invalid α specification")
  @assert(β > 0, "invalid β specification")
  @assert(η > 0, "invalid η specification")
  return function(population::Population, source::Int, target::Int)
    # Ensure source is infectious that target hasn't been previously exposed (SIR model)
    if length(population.events[source][3]) > length(population.events[source][4]) && length(population.events[target][1]) == 0
    # Ensure source is infectious and that target is susceptible (SIS* model)
    #if length(population.events[source][3]) > length(population.events[source][4]) && length(population.events[target][1]) == length(population.events[target][4])
      return α*evaluate(dist, population.history[source][1][1], population.history[target][1][1])^-β
    # Identify an external source and ensure that the target hasn't been previously exposed (SIR model)
    elseif length(population.events[source][1]) > 0 && isnan(population.events[source][1][1]) && length(population.events[target][1]) == 0
    # Identify an external source and ensure that the target is susceptible (SIS* model)
    #elseif isnan(population.events[source][1][1]) && length(population.events[target][1]) == length(population.events[target][4])
      return η
    else
      return 0.
    end
  end
end


"""
Creates generic constant rate function
"""
function create_constantrate(τ::Float64)
  @assert(τ > 0, "invalid τ specification")
  return function(population::Population, individual::Int)
    return τ
  end
end


"""
Generate an array which contains rates (for exponential distribution) for movement from between disease states, and mutation.
`susceptibility_fun` is a function which generates a rate for each pair of target
`substitution_matrix` is a 4x4 array containing single nucleotide polymorphism substitution rates
"""
function create_ratearray(population::Population, susceptibility_fun::Function, substitution_matrix::Array)
  # Set up an array of zeros with rows for each potential source of exposure, for infection, recovery, and mutation at each base location, and columns for each individual...
  rate_array = RateArray(fill(0., (length(population.events)+2+length(population.history[1][2][1]), length(population.events))),
                         fill((0,0,0), (length(population.events)+2+length(population.history[1][2][1]), length(population.events))))

  # Define events
  for r in 1:size(rate_array.events, 1)
    for c in 1:size(rate_array.events, 2)
      if r <= size(rate_array.events, 2)
        # Exposure event (from susceptible to exposed state)
        rate_array.events[r,c] = (1,c,r)
      elseif r == size(rate_array.events, 2)+1
        # Symptom event (from exposed to infectious state)
        rate_array.events[r,c] = (2,c,0)
      elseif r == size(rate_array.events, 2)+2
        # Recovery event (from infectious to recovered or susceptible* state)
        rate_array.events[r,c] = (3,c,0)
      else
        # Mutation event
        rate_array.events[r,c] = (4,c,r-size(rate_array.events, 2)-2)
      end
    end
  end

  # External exposure rate
  for i = 2:size(rate_array.rates,2)
    rate_array.rates[1,i] = susceptibility_fun(population, 1, i)
  end

  # Mutation of external pathogen
  substitution_matrix[[1,6,11,16]] = 0.
  rate_ref = sum(substitution_matrix,2)[:]
  rate_array.rates[length(population.events)+3:end,1] = rate_ref[convert(Vector{Int64}, population.history[1][2][1])]

  # Return RateArray
  return rate_array
end


"""
One event occurs, and appropriate updates are made to the RateArray and Population
"""
function onestep!(rate_array::RateArray, population::Population, susceptibility_fun::Function, latency_fun::Function, recovery_fun::Function, substitution_matrix::Array)
  rate_total = cumsum(rate_array.rates[:])

  if isinf(rate_total[end])
    increment = 0.
    event = rate_array.events[findfirst(isinf(rate_total))]

  else
    increment = rand(Exponential(1/rate_total[end]))
    event = rate_array.events[findfirst(rate_total .> rand()*rate_total[end])]
  end

  # Update population timeline
  push!(population.timeline[1], population.timeline[1][end]+increment)
  push!(population.timeline[2], event)

  if event[1] == 1
    # S => E
    # Update rates - clear all individual specific rates
    rate_array.rates[:,event[2]] = 0.
    # Update rates - latency
    rate_array.rates[size(rate_array.rates,2)+1, event[2]] = latency_fun(population, event[2])
    # Update population - exposure time
    push!(population.events[event[2]][1], population.timeline[1][end])
    # Update population - exposure source
    push!(population.events[event[2]][2], event[3])
    # Update population - sequence
    population.history[event[2]][2] = [population.history[event[3]][2][end]]
    # Update population - sequence time
    push!(population.events[event[2]][6], population.timeline[1][end])
    # Update rates - mutation rates
    substitution_matrix[[1,6,11,16]] = 0.
    rate_ref = sum(substitution_matrix,2)[:]
    rate_array.rates[length(population.events)+3:end, event[2]] = rate_ref[convert(Vector{Int64}, population.history[event[2]][2][1])]

  elseif event[1] == 2
    # E => I
    # Update rates - clear latency
    rate_array.rates[size(rate_array.rates, 2)+1, event[2]] = 0.
    # Update rates - recovery
    rate_array.rates[size(rate_array.rates, 2)+2, event[2]] = recovery_fun(population, event[2])
    # Update population - infection time
    push!(population.events[event[2]][3], population.timeline[1][end])
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
    push!(population.events[event[2]][4], population.timeline[1][end])
    # Update susceptibility* (under SIR framework, susceptibilities of 0 will be generated)
    for i = 1:size(rate_array.rates,2)
      rate_array.rates[i, event[2]] = susceptibility_fun(population, i, event[2])
    end

  else
    # Mutation
    # Update population - sequence
    substitution_matrix[[1,6,11,16]] = 0.
    population.history[event[2]][2] = push!(population.history[event[2]][2], population.history[event[2]][2][end])
    population.history[event[2]][2][end][event[3]] = convert(Nucleotide2bitBase, findfirst(rand(Multinomial(1, substitution_matrix[:,convert(Int64, population.history[event[2]][2][end][event[3]])][:]/sum(substitution_matrix[:,convert(Int64, population.history[event[2]][2][end][event[3]])][:])))))
    # Update population - sequence time
    push!(population.events[event[2]][6], population.timeline[1][end])
    # Update rates - mutation rates
    rate_ref = sum(substitution_matrix,2)[:]
    rate_array.rates[size(rate_array.rates,2)+2+event[3], event[2]] = rate_ref[convert(Int64, population.history[event[2]][2][end][event[3]])]
  end
  return rate_array, population
end


"""
Gather surveillance data on specific individuals in a population, with an exponentially distributed detection lag with rate ν
"""
function surveil(population::Population, ν::Float64)
  exposed_actual = fill(NaN, length(population.events)-1)
  infectious_actual = fill(NaN, length(population.events)-1)
  infectious_observed = fill(NaN, length(population.events)-1)
  removed_actual = fill(NaN, length(population.events)-1)
  removed_observed = fill(NaN, length(population.events)-1)
  covariates_actual = fill(fill(NaN, length(population.history[2][1][1])),  length(population.events)-1)
  covariates_observed = fill(fill(NaN, length(population.history[2][1][1])),  length(population.events)-1)
  seq_actual = convert(Vector{Any}, fill(NaN, length(population.events)-1))
  seq_observed = convert(Vector{Any}, fill(NaN, length(population.events)-1))

  for i = 2:length(population.events)
    # Initial conditions (assumed to be constant and observed without error)
    covariates_actual[i-1] = population.history[i][1][1]
    covariates_observed[i-1] = population.history[i][1][1]

    # Exposure time (unobservable)
    if length(population.events[i][1]) > 0
      exposed_actual[i-1] = population.events[i][1][1]
    end

    # Infectious time (observed with latency)
    if length(population.events[i][3]) > 0
      infectious_actual[i-1] = population.events[i][3][1]
      seq_actual[i-1] = convert(Vector{Int64}, population.history[i][2][find(infectious_actual[i-1] .>= population.events[i][6])[end]])
      if ν < Inf
        infectious_observed[i-1] = infectious_actual[i-1] + rand(Exponential(1/ν))
      elseif ν == Inf
        infectious_observed[i-1] = infectious_actual[i-1]
      end

      if length(population.events[i][4]) > 0 && infectious_observed[i-1] >= population.events[i][4][1]
        infectious_observed[i-1] = NaN
      else
        seq_observed[i-1] = convert(Vector{Int64}, population.history[i][2][find(infectious_observed[i-1] .>= population.events[i][6])[end]])
      end
    end

    # Removal time (observed with latency)
    if length(population.events[i][4]) > 0
      removed_actual[i-1] = population.events[i][4][1]
      if !isnan(infectious_observed[i-1])
        if ν < Inf
          removed_observed[i-1] = removed_actual[i-1] + rand(Exponential(1/ν))
        elseif ν == Inf
          removed_observed[i-1] = removed_actual[i-1]
        end
      end
    end
  end
  if exposed_actual == infectious_actual
    return SIR_actual(infectious_actual, removed_actual, covariates_actual, seq_actual), SIR_observed(infectious_observed, removed_observed, covariates_observed, seq_observed)
  else
    return SEIR_actual(exposed_actual, infectious_actual, removed_actual, covariates_actual, seq_actual), SEIR_observed(infectious_observed, removed_observed, covariates_observed, seq_observed)
  end
end


surveil(population::Population) = surveil(population, Inf)
