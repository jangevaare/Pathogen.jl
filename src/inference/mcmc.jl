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


# function Array(riskparams::RiskParameters)
#   return [riskparams.sparks;
#           riskparams.susceptibility;
#           riskparams.transmissibility;
#           riskparams.infectivity;
#           riskparams.latency;
#           riskparams.removal]
# end


# function Array(trace::PathogenTrace)
#   return [Array(trace.riskparameters[i]) for i = 1:length(trace.riskparameters)]
# end

function length(trace::PathogenTrace)
  return length(trace.riskparameters)
end

function length(riskparameters::RiskParameters)
  return sum([length(riskparameters.sparks);
              length(riskparameters.susceptibility);
              length(riskparameters.transmissibility);
              length(riskparameters.infectivity);
              length(riskparameters.latency);
              length(riskparameters.removal)])
end


function size(trace::PathogenTrace)
  return (length(trace), length(trace[1]))
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


function getindex(x::Array{RiskParameters, 1}, i, j)
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


function Array(trace::PathogenTrace)
  dims = size(trace)
  return [trace.riskparameters[i, j] for i = 1:dim[1], j = 1:dim[2]]
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
Transition kernel
"""
type TransitionKernel
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  latency::Vector{Float64}
  removal::Vector{Float64}

  function(priors::RiskParameterPriors)
    sparks = Float64[]
    susceptibility = Float64[]
    transmissibility = Float64[]
    infectivity = Float64[]
    latency = Float64[]
    removal = Float64[]
    for i in RiskParameterPriors.sparks
      push!(sparks, var(i)*2.38^2)
    end
    for i in RiskParameterPriors.susceptibility
      push!(susceptibility, var(i)*2.38^2)
    end
    for i in RiskParameterPriors.transmissibility
      push!(transmissibility, var(i)*2.38^2)
    end
    for i in RiskParameterPriors.infectivity
      push!(infectivity, var(i)*2.38^2)
    end
    for i in RiskParameterPriors.latency
      push!(latency, var(i)*2.38^2)
    end
    for i in RiskParameterPriors.removal
      push!(removal, var(i)*2.38^2)
    end
    dims = sum([length(sparks);
                length(susceptibility);
                length(transmissibility);
                length(infectivity);
                length(latency);
                length(removal)])

    sparks /= dims
    susceptibility /= dims
    transmissibility /= dims
    infectivity /= dims
    latency /= dims
    removal /= dims

    return new(sparks,
               susceptibility,
               transmissibility,
               infectivity,
               latency,
               removal)
  end

  function(paramstrace::Vector{RiskParameters})

end



"""
Generate proposal
"""
function propose(currentstate::PathogenIteration)
  newstate
  return newstate
end


"""
Phylodynamic ILM MCMC
"""
function mcmc(n::Int64,
              trace::PathogenTrace,
              eventpriors::EventPriors,
              riskparameterpriors::RiskParameterPriors,
              riskfuncs::RiskFunctions,
              population::DataFrame,
              events::Events)
  progressbar = Progress(n, 5, "Performing $n MCMC iterations...", 25)
  for i = 1:n
    next!(progressbar)
    proposal = PathogenProposal


  end
  return trace
end
