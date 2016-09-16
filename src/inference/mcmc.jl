"""
MCMC iteration
"""
type PathogenIteration
  risk_parameters::RiskParameters
  events::Events
  network::Network
end


"""
MCMC trace
"""
type PathogenTrace
  risk_parameters::Vector{RiskParameters}
  events::Vector{Events}
  network::Vector{Network}
end


function push!(trace::PathogenTrace, iteration::PathogenIteration)
  push!(trace.risk_parameters, iteration.risk_parameters)
  push!(trace.events, iteration.events)
  push!(trace.network, iteration.network)
  return trace
end


function append!(trace1::PathogenTrace, trace2::PathogenTrace)
  append!(trace1.risk_parameters, trace2.risk_parameters)
  append!(trace1.events, trace2.events)
  append!(trace1.network, trace2.network)
end


PathogenProposal = PathogenIteration
