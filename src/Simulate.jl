"""
simulate.jl - pathogen evolution and transmission dynamic simulation tools
Justin Angevaare
May 2015
"""

function CreatePopulation(init_seq::Nucleotide2bitSeq, init_var::Array)
"""
Create an infection database.
`init_seq` is assigned to the "external" infection source.
Each row of the `init_array` is assigned to an individual
"""
  # exposure times, exposure source, infection times, recovery times, covariate times, sequence times
  events = Array[Array[[NaN], [NaN],  [NaN],  [NaN],  [NaN], [0]],
                 Array[[],    [],     [],     [],     [0],   []]]

  # covariate history, sequence history
  history = Array[Array[[[fill(NaN, length(init_var[1,:]))]],[init_seq]],
                  Array[[[init_var[1,:]]],[]]]

  # push individuals to these arrays.
  for r = 2:size(init_var,1)
    push!(events, Array[[r], [], [], [], [], [0], []])
    push!(history, Array[[[init_var[r,:]]],[]])
  end

  # save as a population object type
  return population(events, history)
end

function CreatePowerLaw(α::Float64, β::Float64, γ::Float64, η::Float64, dist=Euclidean())
  """
  This function creates a full parameterized power law function
  """
  @assert(α > 0, "invalid α specification")
  @assert(β > 0, "invalid β specification")
  @assert(γ > 0, "invalid γ specification")
  @assert(η > 0, "invalid η specification")

  return function(population::population, source::Vector, target::Vector)
    """
    This simple `SusceptibilityFunction` returns the rate parameter for a `target` individuals from `source` individuals using the power law kernel with parameters α and β. Location must be specified with matching but arbitrary dimensions for each individual; specifically, each individual is represented by a column in an array. Distance by default is Euclidean, but any of the distance calculations in the Distance.jl package may be used.

    A zero distance is assigned a rate of γ, and a NaN distance (external source of infection) is assigned a rate of η.

    It's important to note that this function does not check the disease status of any individuals.

    This function also serves as a model for any user defined
    """
    # `source` individuals as rows and `target` individuals as columns
    rates = α*pairwise(dist, population.history[source][1][end], population.history[target][1][end]).^-β

    # assign rate of γ to for a `target` sharing a location with a `source`
    rates[rates .== 0] = γ

    # assign rate of γ for external sources of infection
    rates[rates .== NaN] = η

    return rates
  end
end

function CreateConstantRate(τ::Float64)
  """
  Creates generic constant rate function
  """
  return function(individuals::Array)
    return fill(τ, shape(individuals, 2))
  end
end

function CreateRateArray(population::population, SusceptibilityFunction::Function, LatencyFunction::Function, RecoveryFunction::Function, SubstitutionMatrix::Array)
  """
  Generate an array which contains rates (for exponential distribution) for movement from between disease states, and mutation.
  `SusceptibilityFunction` is a function which generates a rate for each pair of target
  `SubstitutionMatrix` is a 4x4 array containing single nucleotide polymorphism substitution rates
  """
  # Set up an array of zeros with rows for each potential source of exposure, for infection, recovery, and mutation at each base location, and columns for each individual...
  RateArray = fill(0. (length(population.events)+1+1+length(population.history[1][2][1], length(population.events))))

  # Exposure rate from external source...
  for i = 2:shape(RateArray,2)
    RateArray[1,i] = SusceptibilityFunction(population, 1, i)
  end

  # External source mutation
  RateRef = sum(SubstitutionMatrix,2)
  NucleotideRef = nucleotide2bit("AGCU")
  for i = 1:length(population.history[1][2][1])
    RateArray[length(population.events)+1+1+i,1] = RateRef[population.history[1][2][1][i] .== NucleotideRef]
  end
  return RateArray
end


