function iterate!(mc::MarkovChain{T},
                  mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1) where T <: EpidemicModel
  pmeter = Progress(n, "MCMC progress")
  for i = 1:n
    next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    next!(pmeter)
    @logmsg LogLevel(-5000) "MCMC progress" progress = i/n
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
                  event_batches::Int64=1) where T <: EpidemicModel
  for i = 1:n
    next!(mc, mcmc, Σ, σ, condition_on_network = condition_on_network, event_batches = event_batches)
    put!(progress_channel, true)
  end
  return mc
end

function iterate!(mcmc::MCMC{T},
                  n::Int64,
                  Σ::Array{Float64, 2},
                  σ::Float64;
                  condition_on_network::Bool=false,
                  event_batches::Int64=1) where T <: EpidemicModel
  pmeter = Progress(n*length(mcmc.markov_chains), "MCMC progress")
  pchannel = RemoteChannel(()->Channel{Bool}(10), 1)
  mc_Futures = Future[]
  for mc in mcmc.markov_chains
    @debug "Starting MCMC..."
    push!(mc_Futures, @spawn iterate!(mc, mcmc, n, Σ, σ, pchannel, condition_on_network = condition_on_network, event_batches = event_batches))
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
