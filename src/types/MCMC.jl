mutable struct Trace{T <: EpidemicModel}
  iterations::Int64
  events::Vector{Events{T}}
  network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters}
  log_posterior::Vector{Float64}
end

mutable struct MCMC{T <: EpidemicModel}
  observations::EventObservations{T}
  population::DataFrame
  risk_functions::RiskFunctions{T}
  risk_priors::RiskPriors{T}
  traces::Vector{Trace{T}}

  function MCMC(obs::EventObservations{T},
                pop::DataFrame,
                rf::RiskFunctions{T},
                rp::RiskPriors{T}) where T <: EpidemicModel
    return new{T}(obs, pop, rf, rp, Trace{T}[])
  end
end
