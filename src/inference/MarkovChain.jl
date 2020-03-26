mutable struct MarkovChain{S <: DiseaseStateSequence, M <: ILM}
  iterations::Int64
  events::Vector{Events{M}}
  transmission_network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{M}}
  substitution_model::Union{Nothing, Vector{NucleicAcidSubstitutionModel}}
  log_posterior::Vector{Float64}
  Σrp::Union{Nothing, OnlineStats.CovMatrix}
  Σsm::Union{Nothing, OnlineStats.CovMatrix}

  function MarkovChain(e::Events{S},
                       n::TransmissionNetwork,
                       rp::RiskParameters{S},
                       lp::Float64) where {
                       S <: DiseaseStateSequence,
                       M <: TNILM}
    return new{S, M}(0, [e], [n], [rp], nothing, [lp], OnlineStats.CovMatrix(), nothing)
  end

  function MarkovChain(e::Events{S},
                       n::TransmissionNetwork,
                       rp::RiskParameters{S},
                       sm::NucleicAcidSubstitutionModel,
                       lp::Float64) where {
                       S <: DiseaseStateSequence,
                       M <: PhyloILM}
    return new{S, M}(0, [e], [n], [rp], [sm], [lp], OnlineStats.CovMatrix(), OnlineStats.CovMatrix())
  end
end

function Base.length(x::MarkovChain)
  return x.iterations
end

function Base.show(io::IO, x::MarkovChain{S, M}) where {
                   S <: DiseaseStateSequence,
                   M <: ILM}
  return print(io, "$S $M Markov chain (iterations = $(length(x)))")
end

function TransmissionNetworkDistribution(x::MarkovChain)
  return TNDistribution(x.transmission_network)
end

function TransmissionNetworkDistribution(iter, x::MarkovChain)
  return TNDistribution(x.transmission_network[iter])
end