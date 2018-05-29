mutable struct Trace{T <: EpidemicModel}
  iterations::Int64
  events::Vector{Events{T}}
  network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{T}}
  log_posterior::Vector{Float64}
end
