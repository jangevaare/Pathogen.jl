# type Rates
#   external = Array{Nullable{Float64}, 1}
#   internal = Array{Nullable{Float64}, 2}
#   infection = Array{Nullable{Float64}, 1}
#   removal = Array{Nullable{Float64}, 1}
# end
#
#
# function getindex(x::Rates, i::Real)
#   if i == 1
#     return x.external
#   elseif i == 2
#     return x.internal
#   elseif i == 3
#     return x.infection
#   elseif i == 4
#     return x.removal
#   else
#     throw(BoundsError)
#   end
# end

"""
An array which stores information for simulation purposes
"""
type Rates
  rates::Vector{Array{Float64}}
  mask::Vector{Array{Bool}}

  function Rates(individuals::Int64)
    rates = Array{Float64}[]
    mask = Array{Bool}[]
    # External exposure rates
    push!(rates, fill(0., individuals))
    push!(mask, fill(false, individuals))
    # Internal exposure rates
    push!(rates, fill(0., (individuals, individuals)))
    push!(mask, fill(false, (individuals, individuals)))
    # Infection rates
    push!(rates, fill(0., individuals))
    push!(mask, fill(false, individuals))
    # Removal rates
    push!(rates, fill(0., individuals))
    push!(mask, fill(false, individuals))
    return new(rates, mask)
  end
end

function getindex(x::Rates, i::Real, j::Real)
  return x.mask[i][j] * x.rates[i][j]
end


function getindex(x::Rates, i::Real)
  return x.mask[i] .* x.rates[i]
end
