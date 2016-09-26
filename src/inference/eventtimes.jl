"""
Event time augmentation
"""
function propose(individuals::Vector{Int64},
                 events::Events,
                 eventpriors::EventPriors,
                 network::Network)
  proposal = events
  for i in individuals
    pathfrom = pathwayfrom(i, network)
    pathto = pathwayto(i, network)
    # Exposure time
    if !isnan(events.exposed[i])
      if length(pathto) > 2
        exposure_lb = events.infected[pathto[2]]
        if isnan(events.removed[pathto[2]])
          exposure_ub = events.infected[i]
        else
          exposure_ub = minimum([events.infected[i]; events.removed[pathto[2]]])
        end
      else
        exposure_lb = 0.
        exposure_ub = events.infected[i]
      end
      proposal.exposed[i] = rand(Truncated(eventpriors.exposed[i], exposure_lb, exposure_ub))
    end
    # Infection time
    if !isnan(events.infected[i])
      infection_lb = proposal.exposed[i]
      infection_ub = minimum(proposal.exposed[pathfrom[2:end]], Inf)
      proposal.infected[i] = rand(Truncated(eventpriors.infected[i], infection_lb, infection_ub))
    end
    # Removal time
    if !isnan(events.removed[i])
      removal_lb = maximum(proposal.exposed[pathfrom[2:end]])
      removal_ub = Inf
      proposal.removed[i] = rand(Truncated(eventpriors.removed[i], removal_lb, removal_ub))
    end
  end
  return proposal
end
