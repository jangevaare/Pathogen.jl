function JC69(θ::Tuple, time)
  """
  Returns the JC69 (Jukes and Cantor, 1969) transition matrix with parameter μ, the overall substitution rate. This model assumes equal base frequencies and mutation rates.
  http://en.wikipedia.org/wiki/Models_of_DNA_evolution#Most_common_models_of_DNA_evolution
  """
  μ = θ[1]
  p_0 = 0.25 + 0.75*e^(-time*μ)
  p_1 = 0.25 - 0.25*e^(-time*μ)
  return [[p_0 p_1 p_1 p_1]
          [p_1 p_0 p_1 p_1]
          [p_1 p_1 p_0 p_1]
          [p_1 p_1 p_1 p_0]]
end

function K80(θ::Tuple, time)
  """
  Returns the K80 (Kimura et al., 1980) transition matrix with parameters α and β; the transition and transversion rate parameters respectively. This model assumes equal base frequencies and unique transition and transversion rates
  http://en.wikipedia.org/wiki/Models_of_DNA_evolution#Most_common_models_of_DNA_evolution
  """
  α = θ[1]
  β = θ[2]
  p_0 = 0.25 + 0.25*e^(-4*β*time) + 0.5*e^(-2*(α + β)*time)
  p_1 = 0.25 + 0.25*e^(-4*β*time) - 0.5*e^(-2*(α + β)*time)
  p_2 = 0.25 - 0.25*e^(-4*β*time)
  return [[p_0 p_1 p_2 p_2]
          [p_1 p_0 p_2 p_2]
          [p_2 p_2 p_0 p_1]
          [p_2 p_2 p_1 p_0]]
end

function F81(θ::Tuple, time)
   """
  Returns the F81 (Felsenstein, 1980) transition matrix with parameters μ, π_A, π_G, π_T, and π_C, which are the overall substitution rate and the base specific frequencies respectively.
  http://en.wikipedia.org/wiki/Models_of_DNA_evolution#Most_common_models_of_DNA_evolution
  """
