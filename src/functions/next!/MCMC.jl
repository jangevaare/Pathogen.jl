function next!(mcmc::MCMC{T},
               Σ::Array{Float64, 2},
               σ::Float64) where T <: EpidemicModel
  @simd for mc in mcmc.markov_chains
    next!(mc, mcmc, Σ, σ)
  end
  return mcmc
end
