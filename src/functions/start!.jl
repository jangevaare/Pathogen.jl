function start!(mcmc::MCMC{T},
                markov_chains::Int64;
                attempts::Int64=1000) where T <: EpidemicModel
  mc_Futures = Future[]
  for i = 1:markov_chains
    @debug "creating Markov chain $i initialization Future"
    push!(mc_Futures, @spawn initialize(MarkovChain, mcmc, attempts=attempts))
  end
  @debug "$markov_chains Markov chain initialization Futures successfully created"
  for i = 1:markov_chains
    @debug "Waiting for Markov chain $i initialization to complete"
    wait(mc_Futures[i])
  end
  mcmc.markov_chains = [fetch(mc) for mc in mc_Futures]
  @debug "$markov_chains Markov chains successfully initialized"
  return mcmc
end

function start!(mcmc::MCMC{T};
                attempts::Int64=1000) where T <: EpidemicModel
  push!(mcmc.markov_chains, initialize(MarkovChain, mcmc, attempts=attempts))
  @debug "1 Markov Chain successfully initialized"
  return mcmc
end
