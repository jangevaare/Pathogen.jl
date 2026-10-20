function start!(
  rng::AbstractRNG,
  mcmc::MCMC{S, M},
  n_chains::Int64;
  attempts::Int64=5000) where {
  S <: DiseaseStateSequence,
  M <: ILM}
  if n_chains == 1
    return start!(rng, mcmc; attempts=attempts)
  else
    pmeter = Progress(n_chains * attempts, "Initialization progress")
    pchannel = RemoteChannel(()->Channel{Bool}(10), 1)
    mc_Futures = Future[]
    for i = 1:n_chains
      @logmsg LogLevel(-3000) "Creating Markov chain $i initialization Future"
      chain_rng = MersenneTwister(rand(rng, UInt64))
      push!(mc_Futures, @spawn try
        initialize(chain_rng, MarkovChain, mcmc, pchannel, attempts=attempts)
      finally
        put!(pchannel, false)
      end)
    end
    @logmsg LogLevel(-2000) "$n_chains Markov chain initialization Futures successfully created"
    completed = 0
    while completed < length(mc_Futures)
      if take!(pchannel)
        next!(pmeter)
      else
        completed += 1
      end
    end
    close(pchannel)
    finish!(pmeter)
    markov_chains = [fetch(i) for i in mc_Futures]
    append!(mcmc.markov_chains, markov_chains[nothing.!==markov_chains])
    @logmsg LogLevel(-3000) "$(sum(nothing.!==markov_chains)) Markov chains successfully initialized"
    return mcmc
  end
end

function start!(
  rng::AbstractRNG,
  mcmc::MCMC{S, M};
  attempts::Int64=5000) where {
  S <: DiseaseStateSequence,
  M <: ILM}
  mc = initialize(rng, MarkovChain, mcmc, attempts=attempts)
  if mc !== nothing
    push!(mcmc.markov_chains, mc)
  else
    @error "Failed to initialize a Markov chain"
  end
  return mcmc
end


start!(mcmc::MCMC, args...; kwargs...) = start!(_default_rng(), mcmc, args...; kwargs...)
