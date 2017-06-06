"""
MCMC iteration
"""
type PathogenIteration
  riskparameters::RiskParameters
  substitutionmodel::SubstitutionModel
  events::Events
  network::Network
  logposterior::Float64
end
