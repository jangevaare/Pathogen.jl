function start!(mcmc::MCMC{T},
                markov_chains::Int64;
                attempts::Int64=1000) where T <: EpidemicModel
  mc_Futures = Future[]
  for i = 1:markov_chains
    @logmsg LogLevel(-3000) "Creating Markov chain $i initialization Future"
    push!(mc_Futures, @spawn initialize(MarkovChain, mcmc, attempts=attempts))
  end
  @logmsg LogLevel(-2000) "$markov_chains Markov chain initialization Futures successfully created"
  for i = 1:markov_chains
    @debug "Waiting for Markov chain $i initialization to complete"
    wait(mc_Futures[i])
  end
  mcmc.markov_chains = [fetch(mc) for mc in mc_Futures]
  @logmsg LogLevel(-3000) "$markov_chains Markov chains successfully initialized"
  return mcmc
end

function start!(mcmc::MCMC{T};
                attempts::Int64=1000) where T <: EpidemicModel
  push!(mcmc.markov_chains, initialize(MarkovChain, mcmc, attempts=attempts))
  @debug "1 Markov Chain successfully initialized"
  return mcmc
end
