function start!(mcmc::MCMC{T};
  attempts::Int64=1000,
  markov_chains::Int64=1) where T <: EpidemicModel
  for i = 1:markov_chains
    push!(mcmc.markov_chains, initialize(MarkovChain, mcmc, attempts=attempts))
  end
  return mcmc
end
