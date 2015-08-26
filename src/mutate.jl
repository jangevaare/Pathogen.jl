"""
mutate.jl
"""

function jc69(θ::Tuple)
  """
  Returns the JC69 Q matrix with parameter λ
  """
  λ = θ[1]
  return [[0 λ λ λ]
          [λ 0 λ λ]
          [λ λ 0 λ]
          [λ λ λ 0]]
end

# function k80(θ::Tuple)
#   """
#   Returns the K80 (Kimura et al., 1980) transition matrix with parameters α and β; the transition and transversion rate parameters respectively. This model assumes equal base frequencies and unique transition and transversion rates
#   http://en.wikipedia.org/wiki/Models_of_DNA_evolution#Most_common_models_of_DNA_evolution
#   """
#   α = θ[1]
#   β = θ[2]
# end

# function hky85(θ::Tuple)
#   """
#   Returns the HKY85 (Hasegawa, Kishino and Yano 1985) transition matrix with parameters α, β, π_A, π_T, π_C, and π_G, which are the the transition and transversion rate parameters and the base specific frequencies respectively.
#   """
#   α    = θ[1]
#   β    = θ[2]
#   π_A  = θ[3]
#   π_T  = θ[4]
#   π_C  = θ[5]
#   π_G  = θ[6]
# end
