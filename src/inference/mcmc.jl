"""
MCMC iteration
"""
type PathogenIteration
  riskparameters::RiskParameters
  events::Events
  network::Network
  logposterior::Float64
end


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


PathogenProposal = PathogenIteration


function length(x::PathogenTrace)
  return length(x.riskparameters)
end


function length(x::RiskParameters)
  return sum([length(x.sparks);
              length(x.susceptibility);
              length(x.transmissibility);
              length(x.infectivity);
              length(x.latency);
              length(x.removal)])
end


function size(x::Vector{RiskParameters})
  return (length(x), length(x[1]))
end


function getindex(x::RiskParameters, i)
  inds = cumsum([length(x.sparks);
                 length(x.susceptibility);
                 length(x.transmissibility);
                 length(x.infectivity);
                 length(x.latency);
                 length(x.removal)])
  riskfunc = findfirst(i <= inds)
  if riskfunc == 0
    @error("BoundsError")
  end
  return x.(fieldnames(x)[riskfunc])[end - (inds[riskfunc] - i)]
end


function Vector(x::RiskParameters)
  return [x[i] for i = 1:length(x)]
end


function getindex(x::Vector{RiskParameters}, i, j)
  inds = cumsum([length(x[i].sparks);
                 length(x[i].susceptibility);
                 length(x[i].transmissibility);
                 length(x[i].infectivity);
                 length(x[i].latency);
                 length(x[i].removal)])
  riskfunc = findfirst(j <= inds)
  if riskfunc == 0
    @error("BoundsError")
  end
  return x[i].(fieldnames(x[i])[riskfunc])[end - (inds[riskfunc] - j)]
end


function Array(x::Vector{RiskParameters})
  dims = size(x)
  return [x[i, j] for i = 1:dim[1], j = 1:dim[2]]
end


"""
Metropolis-Hastings rejection using log posteriors
"""
function MHreject(lp1::Float64, lp2::Float64)
  reject = true
  if lp1 >= lp2
    reject = false
  elseif exp(lp1 - lp2) >= rand()
    reject = false
  end
  return reject
end


"""
Metropolis-Hastings acceptance using log posteriors
"""
function MHaccept(lp1::Float64, lp2::Float64)
  return !MHreject(lp1, lp2)
end



"""
Generate the variance-covariance matrix for a MvNormal transition kernel based
upon prior distributions
"""
function transition_kernel_variance(x::RiskParameterPriors)
  diagonal = Float64[]
  for i in x.sparks
    push!(diagonal, var(i)*2.38^2)
  end
  for i in x.susceptibility
    push!(diagonal, var(i)*2.38^2)
  end
  for i in x.transmissibility
    push!(diagonal, var(i)*2.38^2)
  end
  for i in x.infectivity
    push!(diagonal, var(i)*2.38^2)
  end
  for i in x.latency
    push!(diagonal, var(i)*2.38^2)
  end
  for i in x.removal
    push!(diagonal, var(i)*2.38^2)
  end
  diagonal /= length(diagonal)
  return diagm(diagonal)
end


"""
Adapt the variance-covariance matrix for a MvNormal transition kernel
"""
function transition_kernel_variance(x::Vector{RiskParameters})
  return cov(Array(x))*(2.38^2)/length(x[1])
end


"""
Generate proposal
"""
function propose(currentstate::RiskParameters,
                 transition_kernel_variance::Array{Float64, 2})

  newstate = rand(MvNormal(Vector(currentstate), transition_kernel_variance))
  inds = cumsum([length(x[i].sparks);
                 length(x[i].susceptibility);
                 length(x[i].transmissibility);
                 length(x[i].infectivity);
                 length(x[i].latency);
                 length(x[i].removal)])

  return RiskParameters(newstate[1:inds[1]],
                        newstate[inds[1]+1:inds[2]],
                        newstate[inds[2]+1:inds[3]],
                        newstate[inds[3]+1:inds[4]],
                        newstate[inds[4]+1:inds[5]],
                        newstate[inds[5]+1:inds[6]])
end


"""
Phylodynamic ILM MCMC
"""
function mcmc(n::Int64,
              eventpriors::EventPriors,
              riskparameter_priors::RiskParameterPriors,
              riskfuncs::RiskFunctions,
              population::DataFrame,
              tune=1000::Int64)
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  transition_kernel_var = transition_kernel_variance(riskparameterpriors)
  riskparameter_proposal = rand(riskparameter_priors)
  events_proposal = rand(eventpriors)
  lp1, networkrates = loglikelihood(riskparameter_proposal,
                                    events_proposal,
                                    riskfuncs,
                                    population)
  network_proposal = rand(network)
  trace = PathogenTrace([riskparameter_proposal], [events_proposal], [network_proposal])
  for i = 1:n
    next!(progressbar)
    lp2 = lp1
    if mod(tune, i) == 0
      transition_kernel_var = transition_kernel_variance(trace.riskparameters)
    end
    riskparameter_proposal = propose(trace.riskparameters[end], transition_kernel_var)
    events_proposal = propose(events, trace.network[end])
    lp1, networkrates = loglikelihood(riskparameter_proposal,
                                      events_proposal,
                                      riskfuncs,
                                      population)
    if lp1 < Inf
      lp1 += logprior(riskparameter_proposal, riskparameter_priors)
    end





  end
  return trace
end
