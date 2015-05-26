"""
substitute.jl - substitution rate matrices
Justin Angevaare
May 2015
"""

function JC69(θ::Tuple)
  """
  Returns the JC69 Q matrix with parameter λ
  """
  λ = θ[1]
  return [[NaN λ λ λ]
          [λ NaN λ λ]
          [λ λ NaN λ]
          [λ λ λ NaN]]
end

# function JC69(θ::Tuple, time)
#   """
#   Returns the JC69 (Jukes and Cantor, 1969) transition matrix with parameter μ, the overall substitution rate. This model assumes equal base frequencies and mutation rates.
#   http://en.wikipedia.org/wiki/Models_of_DNA_evolution#Most_common_models_of_DNA_evolution
#   """
#   μ = θ[1]
#   p_0 = 0.25 + 0.75*e^(-time*μ)
#   p_1 = 0.25 - 0.25*e^(-time*μ)
#   return [[p_0 p_1 p_1 p_1]
#           [p_1 p_0 p_1 p_1]
#           [p_1 p_1 p_0 p_1]
#           [p_1 p_1 p_1 p_0]]
# end
#
# function K80(θ::Tuple, time)
#   """
#   Returns the K80 (Kimura et al., 1980) transition matrix with parameters α and β; the transition and transversion rate parameters respectively. This model assumes equal base frequencies and unique transition and transversion rates
#   http://en.wikipedia.org/wiki/Models_of_DNA_evolution#Most_common_models_of_DNA_evolution
#   """
#   α = θ[1]
#   β = θ[2]
#   p_0 = 0.25 + 0.25*e^(-4*β*time) + 0.5*e^(-2*(α + β)*time)
#   p_1 = 0.25 + 0.25*e^(-4*β*time) - 0.5*e^(-2*(α + β)*time)
#   p_2 = 0.25 - 0.25*e^(-4*β*time)
#   return [[p_0 p_1 p_2 p_2]
#           [p_1 p_0 p_2 p_2]
#           [p_2 p_2 p_0 p_1]
#           [p_2 p_2 p_1 p_0]]
# end

# function HKY85(θ::Tuple, time)
#   """
#   Returns the HKY85 (Hasegawa, Kishino and Yano 1985) transition matrix with parameters α, β, π_A, π_G, π_C, and π_T, which are the the transition and transversion rate parameters and the base specific frequencies respectively.
#   """
#   α    = θ[1]
#   β    = θ[2]
#   π_A  = θ[3]
#   π_G  = θ[4]
#   π_C  = θ[5]
#   π_T  = θ[6]
# end
