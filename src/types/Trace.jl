mutable struct Trace{T <: EpidemicModel}
  iterations::Int64
  events::Vector{Event{T}}
  network::Vector{TransmissionNetwork}
  risk_parameters::Vector{RiskParameters{T}}
  log_posterior::Vector{Float64}
end
