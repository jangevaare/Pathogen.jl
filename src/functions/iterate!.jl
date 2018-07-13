function iterate!(mc::MarkovChain{T},
                  mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64) where T <: EpidemicModel
  for i = 1:n
    next!(mc, mcmc, Σ, σ)
  end
  return mc
end

function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64) where T <: EpidemicModel
  mc_Futures = Future[]
  for mc in mcmc.markov_chains
     push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σ, σ))
  end
  for mc in mc_Futures
    wait(mc)
  end
  mcmc.markov_chains = [fetch(i) for i in mc_Futures]
  return mcmc
end
