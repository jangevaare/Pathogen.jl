struct RiskPriors{T} <: AbstractRisk{T}
  sparks::Union{Nothing, Vector{UnivariateDistribution}}
  susceptibility::Union{Nothing, Vector{UnivariateDistribution}}
  infectivity::Union{Nothing, Vector{UnivariateDistribution}}
  transmissibility::Union{Nothing, Vector{UnivariateDistribution}}
  latency::Union{Nothing, Vector{UnivariateDistribution}}
  removal::Union{Nothing, Vector{UnivariateDistribution}}

  RiskPriors{SEIR}(pϵ, ps, pi, pκ, pl, pr) = new{SEIR}(pϵ, ps, pi, pκ, pl, pr)
  RiskPriors{SEI}(pϵ, ps, pi, pκ, pl) = new{SEI}(pϵ, ps, pi, pκ, pl, nothing)
  RiskPriors{SIR}(pϵ, ps, pi, pκ, pr) = new{SIR}(pϵ, ps, pi, pκ, nothing, pr)
  RiskPriors{SI}(pϵ, ps, pi, pκ) = new{SI}(pϵ, ps, pi, pκ, nothing, nothing)
end


function Base.show(io::IO, x::RiskPriors{T}) where T <: EpidemicModel
  return print(io, "$T model risk function priors")
end
