function start!(mcmc::MCMC{T},
                markov_chains::Int64;
                attempts::Int64=1000) where T <: EpidemicModel
  pmeter = Progress(markov_chains*attempts, "Initialization progress")
  pchannel = RemoteChannel(()->Channel{Bool}(10), 1)
  mc_Futures = Future[]
  for i = 1:markov_chains
    @logmsg LogLevel(-3000) "Creating Markov chain $i initialization Future"
    push!(mc_Futures, @spawn initialize(MarkovChain, mcmc, pchannel, attempts=attempts))
  end
  @logmsg LogLevel(-2000) "$markov_chains Markov chain initialization Futures successfully created"
  while !all(isready.(mc_Futures))
    take!(pchannel) && next!(pmeter)
  end
  close(pchannel)
  finish!(pmeter)
  mcmc.markov_chains = [fetch(i) for i in mc_Futures]
  @logmsg LogLevel(-3000) "$markov_chains Markov chains successfully initialized"
  return mcmc
end

function start!(mcmc::MCMC{T};
                attempts::Int64=1000) where T <: EpidemicModel
  return start!(mcmc, 1, attempts)
end
