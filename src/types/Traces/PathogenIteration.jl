"""
MCMC iteration
"""
type PathogenIteration
  riskparameters::RiskParameters
  events::Events
  network::Network
  logposterior::Float64
end
