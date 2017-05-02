"""
MCMC trace
"""
type PathogenTrace
  riskparameters::Vector{RiskParameters}
  events::Vector{Events}
  network::Vector{Network}
  logposterior::Vector{Float64}
end


function push!(trace::PathogenTrace, iteration::PathogenIteration)
  push!(trace.riskparameters, iteration.riskparameters)
  push!(trace.events, iteration.events)
  push!(trace.network, iteration.network)
  push!(trace.logposterior, iteration.logposterior)
  return trace
end


function append!(trace1::PathogenTrace, trace2::PathogenTrace)
  append!(trace1.riskparameters, trace2.riskparameters)
  append!(trace1.events, trace2.events)
  append!(trace1.network, trace2.network)
  append!(trace1.logposterior, trace2.logposterior)
  return trace1
end


function deleteat!(trace::PathogenTrace, inds)
  deleteat!(trace.riskparameters, inds)
  deleteat!(trace.events, inds)
  deleteat!(trace.network, inds)
  deleteat!(trace.logposterior, inds)
  return trace
end


function length(x::PathogenTrace)
  return length(x.riskparameters)
end


function getindex(trace::PathogenTrace, i)
  return PathogenTrace(trace.riskparameters[i],
                       trace.events[i],
                       trace.network[i],
                       trace.logposterior[i])
end


function show(io::IO, object::PathogenTrace)
  print(io, "PathogenTrace object (MCMC iterations: $(length(object)))")
end
