"""
simulate.jl - pathogen evolution and transmission dynamic simulation tools
Justin Angevaare
May 2015
"""

function create_population(init_seq::Nucleotide2bitSeq, init_var::Array)
"""
Create an infection database.
`init_seq` is assigned to the "external" infection source.
Each row of the `init_array` is assigned to an individual
"""
  # identifier, exposure times, exposure source, infection times, recovery times, covariate times, sequence times
  events = Array[Array[[0], [NaN], [NaN],  [NaN],  [NaN],  [NaN], [0]],
                 Array[[1], [],    [],     [],     [],     [0],   []]]

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

#function CreateRateArray(population::population, ExternalPressure, SusceptibilityFunction::Function, LatentPeriod::Float64, InfectiousPeriod::Float64, SubstitutionMatrix::Array)
function CreateRateArray(population::population, SusceptibilityFunction::Function, LatencyFunction::Function, RecoveryFunction::Function, SubstitutionMatrix::Array)
  """
  Generate an array which contains rates (for exponential distribution) for movement from between disease states, and mutation.
  `external_pressure` is a disease exposure rate from outside the described population
  `susceptibility_function` is a user defined function which generates a rate for each pair of individuals based on their covariates and infection history
  `latent_period` is the average length of time spent by an individual in an exposed state (set to 0 if exposed state does not occur)
  `infectious_period` is the average length of time spent by an individual in an infected state (set to 0 if infectious state does not occur)
  `substitution_matrix` is a 4x4 array containing single nucleotide polymorphism substitution rates
  """
  # Set up an array of zeros with rows for each potential source of exposure, for infection, recovery, and mutation at each base location, and columns for each individual...
  RateArray = fill(0. (length(population.events)+1+1+length(population.history[1][2][1], length(population.events))))

  # Exposure rate from external source...
  RateArray[1,2:end] = ExternalPressure

  # External source mutation
  RateRef = sum(SubstitutionMatrix,1)
  NucleotideRef = nucleotide2bit("AGCU")
  for i = 1:length(population.history[1][2][1])
    RateArray[length(population.events)+1+1+i,1] = RateRef[population.history[1][2][1][i] .== NucleotideRef]
  end
  return RateArray
end
