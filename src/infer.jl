"""
Prior distributions vectors for the `RiskParameters`
"""
type RiskPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
  detection::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}
end


"""
MCMC trace
"""
type PathogenTrace
  risk_parameters::Vector{RiskParameters}
  events::Vector{Events}
end


"""
MCMC iteration
"""
type PathogenIteration
  risk_parameters::RiskParameters
  events::Events
end


PathogenProposal = PathogenIteration


function push!(trace::PathogenTrace, iteration::PathogenIteration)
  push!(trace.risk_parameters, iteration.risk_parameters)
  push!(trace.events, iteration.events)
  return trace
end


function append!(trace1::PathogenTrace, trace2::PathogenTrace)
  append!(trace1.risk_parameters, trace2.risk_parameters)
  append!(trace1.events, trace2.events)
end
