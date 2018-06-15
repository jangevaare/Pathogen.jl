mutable struct MarkovChain{T <: EpidemicModel}
  iterations::Int64
  events::Vector{Events{T}}
  network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{T}}
  log_posterior::Vector{Float64}

  function MarkovChain(iterations::Int64,
                       events::Vector{Events{T}},
                       network::Vector{TransmissionNetwork},
                       risk_parameters::Vector{RiskParameters{T}},
                       log_posterior::Vector{Float64}) where T <: EpidemicModel
    return new{T}(iterations, events, network, risk_parameters, log_posterior)
  end

  function MarkovChain{T}() where T <: EpidemicModel
    return new{T}()
  end
end
