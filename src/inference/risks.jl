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
function logprior(riskparams::RiskParameters,
                  riskpriors::RiskParameterPriors)
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
  return diagm(diagonal)
end


"""
Adapt the variance-covariance matrix for a MvNormal transition kernel for
`RiskParameters`
"""
function transition_kernel_variance(x::Vector{RiskParameters})
  return cov(Array(x))*(2.38^2)/length(x[1])
end


"""
Generate a `RiskParameters` proposal using the multivariate normal distribution
as the transition kernel, with a previous set of `RiskParameters` as the mean
vector and a transition kernel variance as the variance-covariance matrix
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
