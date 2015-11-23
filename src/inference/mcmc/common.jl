"""
Calculate the log prior from prior distributions and specified parameter values
"""
function logprior(priors::Priors, params::Vector{Float64}, debug=false::Bool)
  @assert(length(params) == length(fieldnames(priors)),
          "Mismatch between parameter vector and prior")
  lprior = 0.
  for i = 1:length(params)
    lprior += logpdf(priors.(fieldnames(priors)[i]), params[i])
  end
  debug && println("$(typeof(priors)) log prior: $(round(lprior,3))")
  return lprior
end


"""
Randomly generate a parameter vector from specified priors
"""
function randprior(priors::Priors)
  params = [rand(priors.(fieldnames(priors)[1]))]
  for i = 2:length(fieldnames(priors))
    push!(params, rand(priors.(fieldnames(priors)[i])))
  end
  return params
end


"""
A simple function for Metropolis-Hastings rejection using log posteriors
"""
function MHreject(lp1::Float64, lp2::Float64, debug=false::Bool)
  @assert(lp1 < Inf && lp2 < Inf, "Infinite log posterior detected")
  reject = true
  if lp1 >= lp2
    debug && println("MCMC proposal accepted")
    reject = false
  elseif exp(lp1 - lp2) >= rand()
    debug && println("MCMC proposal probabilistically accepted")
    reject = false
  end
  debug && reject && println("MCMC Proposal rejected")
  return reject
end
