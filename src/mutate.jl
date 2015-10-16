"""
Returns the JC69 Q matrix with parameter λ
"""
function jc69(θ::Vector{Float64})
  λ = θ[1]
  return [[-3λ λ λ λ]
          [λ -3λ λ λ]
          [λ λ -3λ λ]
          [λ λ λ -3λ]]
end


"""
Returns the JC69 Q matrix with parameter λ
"""
function jc69q(θ::Vector{Float64})
  λ = θ[1]
  return [[-3λ λ λ λ]
          [λ -3λ λ λ]
          [λ λ -3λ λ]
          [λ λ λ -3λ]]
end


"""
Returns the JC69 P matrix with parameter λ, for time t
"""
function jc69p(θ::Vector{Float64}, t::Float64)
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


#   """
#   Returns the K80 (Kimura et al., 1980) Q matrix with parameters α and β
#   """
# function k80q(θ::Vector{Float64})
#   α = θ[1]
#   β = θ[2]
# end


#   """
#   Returns the K80 (Kimura et al., 1980) P matrix with parameters α and β, for time t
#   """
# function k80p(θ::Vector{Float64}, t::Float64)
#   α = θ[1]
#   β = θ[2]
# end


# function k80p(θ::Vector{Float64})
#   return function(t)
#     k80p(θ, t)
#   end
# end


#   """
#   Returns the HKY85 (Hasegawa, Kishino and Yano 1985) Q matrix with parameters α, β, π_A, π_T, π_C, and π_G
#   """
# function hky85q(θ::Vector{Float64})
#   α    = θ[1]
#   β    = θ[2]
#   π_A  = θ[3]
#   π_T  = θ[4]
#   π_C  = θ[5]
#   π_G  = θ[6]
# end


#   """
#   Returns the HKY85 (Hasegawa, Kishino and Yano 1985) P matrix with parameters α, β, π_A, π_T, π_C, and π_G, for time t
#   """
# function hky85p(θ::Vector{Float64}, t::Float64)
#   α    = θ[1]
#   β    = θ[2]
#   π_A  = θ[3]
#   π_T  = θ[4]
#   π_C  = θ[5]
#   π_G  = θ[6]
# end


# function hky85p(θ::Vector{Float64})
#   return function(t)
#     hky85p(θ, t)
#   end
# end
