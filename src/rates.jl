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
    # Detection rates
    push!(rates, fill(0., individuals))
    push!(mask, fill(false, individuals))
    # Removal rates
    push!(rates, fill(0., individuals))
    push!(mask, fill(false, individuals))
    return new(rates, mask)
  end
end


import Base.getindex


function getindex(x::Rates, i, j)
  return x.mask[i][j] * x.rates[i][j]
end


function getindex(x::Rates, i)
  return x.mask[i] .* x.rates[i]
end
