function iterate!(mc::MarkovChain{T},
                  mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64) where T <: EpidemicModel
  for i = 1:n
    next!(mc, mcmc, Σ, σ)
    @logmsg LogLevel(-5000) "MCMC progress" maxlog=1 core = myid() progress = i/n
  end
  return mc
end

function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64) where T <: EpidemicModel
  mc_Futures = Future[]
  for mc in mcmc.markov_chains
    @debug "Beginning Markov chain iteration"
    push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σ, σ))
  end
  for mc in mc_Futures
    @debug "Waiting for Markov chain iteration to complete"
    wait(mc)
  end
  mcmc.markov_chains = [fetch(i) for i in mc_Futures]
  @info "$n iterations for $(length(mcmc.markov_chains)) Markov chains successfully completed"
  return mcmc
end
