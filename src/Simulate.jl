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
  history = Array[Array[[[fill(NaN, length(init_var[1,:]))]],     [init_seq]],
                  Array[[[init_var[1,:]]],[]]]

  # push individuals to these arrays.
  for r = 2:shape(init_var,1)
    push!(events, Array[[r], [], [], [], [], [0], []])
    push!(history, Array[[[init_var[r,:]]],[]])
  end

  # save as a population object type
  return population(events, history)
end

function rate_array(population::population, external_pressure, susceptibility_function::Function, latent_period::Float64, infectious_period::Float64, substitution_matrix::Array)
  """
  Generate an array which contains rates (for exponential distribution) for movement from between disease states, and mutation.
  `external_pressure` is a disease exposure rate from outside the described population
  `susceptibility_function` is a user defined function which generates a rate for each pair of individuals based on their covariates and infection history
  `latent_period` is the average length of time spent by an individual in an exposed state (set to 0 if exposed state does not occur)
  `infectious_period` is the average length of time spent by an individual in an infected state (set to 0 if infectious state does not occur)
  `substitution_matrix` is a 4x4 array containing single nucleotide polymorphism substitution rates
  """
end

