function simulate!(sim::Simulation{T}; tmax::Float64=Inf, nmax::Int64=1000, pmax::Float64=Inf) where T <: EpidemicModel
  if tmax <= 0.0
    error("The simulation time maximum must be > 0.0")
  elseif nmax < 1
    error("The iteration maximum must be > 0")
  elseif pmax <= 0.0
    error("The processing time maximum must be > 0.0")
  else
    println("Simulating until: \n* a simulation time maximum of $tmax,\n* an iteration maximum of $nmax, or,\n* a processing time maximum of $pmax seconds.\n")
  end
  stoptime = time() + pmax
  while true
    next!(sim)
    if sim.time >= tmax
      println("Simulation completed: time maximum reached")
      break
    elseif sim.iterations >= nmax
      println("Simulation completed: iteration maximum reached")
      break
    elseif pmax < Inf && time() >= stoptime
      println("Simulation stopped: processing time maximum reached")
      break
    end
  end
  return sim
end
