function iterate!(mc::MarkovChain{T},
                  mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  condition_on_network::Bool=false,
                  batch_size::Int64=1,
                  adapt_cov::Int64=100) where T <: EpidemicModel
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  pmeter = Progress(n, "MCMC progress")
  use_adapted_cov = false
  local adapted_cov
  if (adapt_cov != 0 & mc.iterations > adapt_cov) && LinearAlgebra.isposdef(value(mc.cov))
    use_adapted_cov = true
    adapted_cov = OnlineStats.value(mc.cov) * 2.38^2 / length(mc.risk_parameters[1])
  end
  for i = 1:n
    if use_adapted_cov
      next!(mc, mcmc, adapted_cov, σ, condition_on_network = condition_on_network, batch_size = batch_size)
    else
      next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, batch_size = batch_size)
    end
    next!(pmeter)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
    if adapt_cov != 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.cov, convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end]))
      if LinearAlgebra.isposdef(value(mc.cov))
        if !use_adapted_cov
          use_adapted_cov = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_cov = OnlineStats.value(mc.cov) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_cov
      end
    end
  end
  return mc
end

function iterate!(mc::MarkovChain{T},
                  mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64,
                  progress_channel::RemoteChannel;
                  condition_on_network::Bool=false,
                  batch_size::Int64=1,
                  adapt_cov::Int64=100) where T <: EpidemicModel
  use_adapted_cov = false
  local adapted_cov
  if (adapt_cov != 0 & mc.cov.n >= adapt_cov) && LinearAlgebra.isposdef(value(mc.cov))
    use_adapted_cov = true
    adapted_cov = OnlineStats.value(mc.cov) * 2.38^2 / length(mc.risk_parameters[1])
  end
  for i = 1:n
    if use_adapted_cov
      next!(mc, mcmc, adapted_cov, σ, condition_on_network = condition_on_network, batch_size = batch_size)
    else
      next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, batch_size = batch_size)
    end
    put!(progress_channel, true)
    if adapt_cov != 0 && mod(i, adapt_cov) == 0
      OnlineStats.fit!(mc.cov, convert(Array{Float64, 2}, mc.risk_parameters[end-adapt_cov:end]))
      if LinearAlgebra.isposdef(value(mc.cov))
        if !use_adapted_cov
          use_adapted_cov = true
          @debug "Now using adapted covariance matrix for adaptive MCMC on core $(Distributed.myid()) (updated every $adapt_cov iterations)"
        end
        adapted_cov = OnlineStats.value(mc.cov) * 2.38^2 / length(mc.risk_parameters[1])
        @debug "Covariance matrix updated for Adaptive Metropolis-Hastings MCMC" Covariance = adapted_cov
      end
    end
  end
  return mc
end

function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  condition_on_network::Bool=false,
                  batch_size::Int64=1,
                  adapt_cov::Int64=100) where T <: EpidemicModel
  if adapt_cov < 0
    @warn "`adapt_cov` argument indicates the increment in iterations in which the covariance matrix is updated and must be ≧ 0. Setting to 0 for non-adaptive Metropolis-Hastings."
    adapt_cov = 0
  end
  pmeter = Progress(n*length(mcmc.markov_chains), "MCMC progress")
  pchannel = RemoteChannel(()->Channel{Bool}(10), 1)
  mc_Futures = Future[]
  for mc in mcmc.markov_chains
    @debug "Starting MCMC..."
    push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σ, σ, pchannel, condition_on_network = condition_on_network, batch_size = batch_size, adapt_cov = adapt_cov))
  end
  @debug "MCMC in progress..."
  while !all(isready.(mc_Futures))
    take!(pchannel) && next!(pmeter)
  end
  close(pchannel)
  finish!(pmeter)
  mcmc.markov_chains = [fetch(i) for i in mc_Futures]
  @debug "MCMC complete..."
  return mcmc
end

function iterate!(mc::MarkovChain{T},
                  mcmc::MCMC{T},
                  n::Int64,
                  σ::Float64;
                  condition_on_network::Bool=false,
                  batch_size::Int64=1,
                  adapt_cov::Int64=100) where T <: EpidemicModel
  return iterate!(mc,
                  mcmc,
                  n,
                  Matrix(LinearAlgebra.Diagonal([var(mcmc.risk_priors[i]) for i in 1:length(mcmc.risk_priors)])),
                  σ,
                  condition_on_network = condition_on_network,
                  batch_size = batch_size,
                  adapt_cov = adapt_cov)
end

function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  σ::Float64;
                  condition_on_network::Bool=false,
                  batch_size::Int64=1,
                  adapt_cov::Int64=100) where T <: EpidemicModel
  return iterate!(mcmc,
                  n,
                  Matrix(LinearAlgebra.Diagonal([var(mcmc.risk_priors[i]) for i in 1:length(mcmc.risk_priors)]))/(10*length(mcmc.risk_priors)),
                  σ,
                  condition_on_network = condition_on_network,
                  batch_size = batch_size,
                  adapt_cov = adapt_cov)
end
