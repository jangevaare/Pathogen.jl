function start!(mcmc::MCMC{T},
                markov_chains::Int64;
                attempts::Int64=1000,
                debug_level::Int64=0) where T <: EpidemicModel
  mc_Futures = Future[]
  for i = 1:markov_chains
      if debug_level >= 3
        println("start!: creating Markov chain $i initialization Future")
      end
     push!(mc_Futures, @spawn initialize(MarkovChain, mcmc, attempts=attempts))
  end
  if debug_level >= 2
    println("start!: $markov_chains Markov chain initialization Futures successfully created")
  end
  for i = 1:markov_chains
    if debug_level >= 3
      println("start!: waiting for Markov chain $i initialization to complete")
    end
    wait(mc_Futures[i])
  end
  mcmc.markov_chains = [fetch(mc) for mc in mc_Futures]
  if debug_level >= 1
    println("start!: $markov_chains Markov chains successfully initialized")
  end
  return mcmc
end

function start!(mcmc::MCMC{T};
                attempts::Int64=1000,
                debug_level::Int64=0) where T <: EpidemicModel
  push!(mcmc.markov_chains, initialize(MarkovChain, mcmc, attempts=attempts))
  if debug_level >= 1
    println("start!: 1 Markov Chain successfully initialized")
  end
  return mcmc
end
