mutable struct MarkovChain{T <: EpidemicModel}
  iterations::Int64
  events::Vector{Events{T}}
  transmission_network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{T}}
  log_posterior::Vector{Float64}
  cov::OnlineStats.CovMatrix

  function MarkovChain(events::Events{T}, network::TransmissionNetwork, risk_parameters::RiskParameters{T}, lp::Float64) where  T <: EpidemicModel
    return new{T}(0, [events], [network], [risk_parameters], [lp], OnlineStats.CovMatrix())
  end

  function MarkovChain{T}() where  T <: EpidemicModel
    return new{T}(0, Events{T}[], TransmissionNetwork[], RiskParameters{T}[], Float64[], OnlineStats.CovMatrix())
  end
end

function Base.length(x::MarkovChain{T}) where T <: EpidemicModel
  return x.iterations
end

function Base.show(io::IO, x::MarkovChain{T}) where T <: EpidemicModel
  return print(io, "$T model Markov chain (iterations = $(length(x)))")
end
