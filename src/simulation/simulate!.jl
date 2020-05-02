function simulate!(sim::Simulation; tmax::Float64=Inf, nmax::Int64=100000, pmax::Float64=Inf)
  if tmax <= 0.0
    throw(ErrorException("The simulation time maximum must be > 0.0"))
  elseif nmax < 1
    throw(ErrorException("The iteration maximum must be > 0"))
  elseif pmax <= 0.0
    throw(ErrorException("The processing time maximum must be > 0.0"))
  else
    @info "Simulating until: \n* a simulation time maximum of $tmax,\n* an iteration maximum of $nmax, or,\n* a processing time maximum of $pmax seconds.\n"
  end
  stoptime = time() + pmax
  while true
    event = generate(AbstractEvent, sim.event_rates, sim.time)
    if _time(event) > tmax
      sim.time = tmax
      @info "Simulation stopped: simulation time maximum reached"
      break
    else
      update!(sim, event)
    end
    if sim.iterations >= nmax
      @info "Simulation stopped: iteration maximum reached"
      break
    elseif pmax < Inf && time() >= stoptime
      @info "Simulation stopped: processing time maximum reached"
      break
    end
  end
  return sim
end
