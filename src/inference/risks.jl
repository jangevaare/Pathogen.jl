"""
Prior distributions vectors for the `RiskParameters`
"""
type RiskParameterPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}
end


"""
Randomly generate a set of `RiskParameters` from their prior distributions
"""
function rand(riskpriors::RiskParameterPriors)
  sparks = Float64[]
  susceptibility = Float64[]
  transmissibility = Float64[]
  infectivity = Float64[]
  latency = Float64[]
  removal = Float64[]

  for i = 1:length(riskpriors.sparks)
    push!(sparks, rand(riskpriors.sparks[i]))
  end

  for i = 1:length(riskpriors.susceptibility)
    push!(susceptibility, rand(riskpriors.susceptibility[i]))
  end

  for i = 1:length(riskpriors.transmissibility)
    push!(transmissibility, rand(riskpriors.transmissibility[i]))
  end

  for i = 1:length(riskpriors.infectivity)
    push!(infectivity, rand(riskpriors.infectivity[i]))
  end

  for i = 1:length(riskpriors.latency)
    push!(latency, rand(riskpriors.latency[i]))
  end

  for i = 1:length(riskpriors.removal)
    push!(removal, rand(riskpriors.removal[i]))
  end

  return RiskParameters(sparks,
                        susceptibility,
                        transmissibility,
                        infectivity,
                        latency,
                        removal)
end


"""
Calculate the log prior of a set of `RiskParameters`
"""
function logprior(riskpriors::RiskParameterPriors,
                  riskparams::RiskParameters)
  lp = 0.
  for i = 1:length(riskparams.sparks)
    lp += loglikelihood(riskpriors.sparks[i], [riskparams.sparks[i]])
  end
  for i = 1:length(riskparams.susceptibility)
    lp += loglikelihood(riskpriors.susceptibility[i], [riskparams.susceptibility[i]])
  end
  for i = 1:length(riskparams.transmissibility)
    lp += loglikelihood(riskpriors.transmissibility[i], [riskparams.transmissibility[i]])
  end
  for i = 1:length(riskparams.infectivity)
    lp += loglikelihood(riskpriors.infectivity[i], [riskparams.infectivity[i]])
  end
  for i = 1:length(riskparams.latency)
    lp += loglikelihood(riskpriors.latency[i], [riskparams.latency[i]])
  end
  for i = 1:length(riskparams.removal)
    lp += loglikelihood(riskpriors.removal[i], [riskparams.removal[i]])
  end
  return lp
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
  return diagonal
end


"""
Adapt the variance-covariance matrix for a MvNormal transition kernel for
`RiskParameters`
"""
function transition_kernel_variance(x::Vector{RiskParameters})
  kernel_var = cov(Array(x))*(2.38^2)/length(x[1])
  return diag(kernel_var)
end


"""
Generate a `RiskParameters` proposal using the multivariate normal distribution
as the transition kernel, with a previous set of `RiskParameters` as the mean
vector and a transition kernel variance as the variance-covariance matrix
"""
function propose(currentstate::RiskParameters,
                 riskpriors::RiskParameterPriors,
                 variance::Vector{Float64})
                 newstate = currentstate
                 variance_index = 1
  for i = 1:length(riskpriors.sparks)
    lb = support(riskpriors.sparks[i]).lb
    ub = support(riskpriors.sparks[i]).ub
    newstate.sparks[i] = rand(Truncated(Normal(currentstate.sparks[i], variance[variance_index]), lb, ub))
    variance_index += 1
  end
  for i = 1:length(riskpriors.susceptibility)
    lb = support(riskpriors.susceptibility[i]).lb
    ub = support(riskpriors.susceptibility[i]).ub
    newstate.susceptibility[i] = rand(Truncated(Normal(currentstate.susceptibility[i], variance[variance_index]), lb, ub))
    variance_index += 1
  end
  for i = 1:length(riskpriors.transmissibility)
    lb = support(riskpriors.transmissibility[i]).lb
    ub = support(riskpriors.transmissibility[i]).ub
    newstate.transmissibility[i] = rand(Truncated(Normal(currentstate.transmissibility[i], variance[variance_index]), lb, ub))
    variance_index += 1
  end
  for i = 1:length(riskpriors.infectivity)
    lb = support(riskpriors.infectivity[i]).lb
    ub = support(riskpriors.infectivity[i]).ub
    newstate.infectivity[i] = rand(Truncated(Normal(currentstate.infectivity[i], variance[variance_index]), lb, ub))
    variance_index += 1
  end
  for i = 1:length(riskpriors.latency)
    lb = support(riskpriors.latency[i]).lb
    ub = support(riskpriors.latency[i]).ub
    newstate.latency[i] = rand(Truncated(Normal(currentstate.latency[i], variance[variance_index]), lb, ub))
    variance_index += 1
  end
  for i = 1:length(riskpriors.removal)
    lb = support(riskpriors.removal[i]).lb
    ub = support(riskpriors.removal[i]).ub
    newstate.removal[i] = rand(Truncated(Normal(currentstate.removal[i], variance[variance_index]), lb, ub))
    variance_index += 1
  end
  return newstate
end
