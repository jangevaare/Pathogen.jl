mutable struct MarkovChain{T <: EpidemicModel}
  iterations::Int64
  events::Vector{Events{T}}
  network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{T}}
  log_posterior::Vector{Float64}
end

function Base.length(x::MarkovChain{T}) where T <: EpidemicModel
  return x.iterations
end

function Base.show(io::IO, x::MarkovChain{T}) where T <: EpidemicModel
  return print(io, "$T model Markov chain (iterations = $(length(x))")
end
