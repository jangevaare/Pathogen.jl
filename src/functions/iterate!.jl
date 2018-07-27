function iterate!(mc::MarkovChain{T},
                  mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  debug_level::Int64=0) where T <: EpidemicModel
  pm = Progress(n, "Performing $n MCMC iterations (pid = $(myid()))")
  for i = 1:n
    next!(pm)
    next!(mc, mcmc, Σ, σ, debug_level=debug_level)
  end
  return mc
end

function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  debug_level::Int64=0) where T <: EpidemicModel
  mc_Futures = Future[]
  for mc in mcmc.markov_chains
      if debug_level >= 3
        println("iterate!: beginning iteration of Markov chain $i")
      end
     push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σ, σ))
  end
  for mc in mc_Futures
    if debug_level >= 3
      println("iterate!: waiting for iteration of Markov chain $i to complete")
    end
    wait(mc)
  end
  mcmc.markov_chains = [fetch(i) for i in mc_Futures]
  if debug_level >= 1
    println("iterate!: $n iterations for $(length(mcmc.markov_chains)) Markov chains successfully completed")
  end
  return mcmc
end
