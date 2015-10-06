"""
mutate.jl
"""

function jc69(θ::Vector{Float64})
  """
  Returns the JC69 Q matrix with parameter λ
  """
  λ = θ[1]
  return [[-3λ λ λ λ]
          [λ -3λ λ λ]
          [λ λ -3λ λ]
          [λ λ λ -3λ]]
end

function jc69q(θ::Vector{Float64})
  """
  Returns the JC69 Q matrix with parameter λ
  """
  λ = θ[1]
  return [[-3λ λ λ λ]
          [λ -3λ λ λ]
          [λ λ -3λ λ]
          [λ λ λ -3λ]]
end

function jc69p(θ::Vector{Float64}, t::Float64)
  """
  Returns the JC69 P matrix with parameter λ, for time t

  Molecular Evolution: a statistical approach, Z. Yang
  """
  λ = θ[1]
  p0 = 0.25 + 0.75*exp(-t*λ)
  p1 = 0.25 - 0.25*exp(-t*λ)
  return [[p0 p1 p1 p1]
          [p1 p0 p1 p1]
          [p1 p1 p0 p1]
          [p1 p1 p1 p0]]
end

function jc69p(θ::Vector{Float64})
  return function(t)
    jc69p(θ, t)
  end
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
