"""
Prior distributions vectors for the `RiskParameters`
"""
type RiskPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{Float64}
  detection::Vector{Float64}
  removal::Vector{Float64}
end


"""
MCMC trace
"""
type Trace
  risk_parameters::Vector{RiskParameters}
  events::Vector{Events}
  tree::Vector{Tree}
end


"""
MCMC iteration
"""
type Iteration
  risk_parameters::RiskParameters
  events::Events
  tree::Tree
end


Proposal = Iteration


function push!(trace::Trace, iteration::Iteration)
  push!(trace.risk_parameters, iteration.risk_parameters)
  push!(trace.events, iteration.events)
  push!(trace.tree, iteration.tree)
  return trace
end
