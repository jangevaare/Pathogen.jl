"""
propose(network_rates::NetworkRates)

Probablistically generate a network object based on exposure network rates
"""
function propose(network_rates::NetworkRates)
  external_rates = network_rates.external
  internal_rates = network_rates.internal
  external_network = fill(false, length(external_rates))
  internal_network = fill(false, size(internal_rates))
  if !(length(external_rates) == size(internal_rates, 1) == size(internal_rates, 2))
    throw(BoundsError)
  end
  for i = 1:length(external_rates)
    external_total = external_rates[i]
    internal_total = sum(internal_rates[:, i])
    if sum(external_total + internal_total) > 0.
      if rand() < external_total/(external_total + internal_total)
        external_network[i] = true
      else
        source = findfirst(rand(Multinomial(1, internal_rates[:, i]/internal_total)))
        internal_network[source, i] = true
      end
    end
  end
  return Network(external_network, internal_network)
end


"""
propose(individuals::Vector{Int64},
        network::Network,
        network_rates::NetworkRates)

Propose an exposure network based on a previous exposure network and exposure
network rates
"""
function propose(individuals::Vector{Int64},
                 network::Network,
                 network_rates::NetworkRates)
  external = copy(network.external)
  internal = copy(network.internal)
  external_rates = network_rates.external
  internal_rates = network_rates.internal
  for i in individuals
    external[i] = false
    internal[:, i] = false
    external_total = external_rates[i]
    internal_total = sum(internal_rates[:, i])
    if sum(external_total + internal_total) > 0.
      if rand() < external_total/(external_total + internal_total)
        external[i] = true
      else
        source = findfirst(rand(Multinomial(1, internal_rates[:, i]/internal_total)))
        internal[source, i] = true
      end
    end
  end
  return Network(external, internal)
end


"""
propose(i::Int64,
        network::Network,
        network_rates::NetworkRates)

Propose an exposure network based on a previous exposure network and exposure
network rates
"""
function propose(i::Int64,
                 network::Network,
                 network_rates::NetworkRates)
  external = network.external
  internal = network.internal
  external_rates = network_rates.external
  internal_rates = network_rates.internal
  external_total = external_rates[i]
  internal_total = sum(internal_rates[:, i])
  if sum(external_total + internal_total) > 0.
    external[i] = false
    internal[:, i] = false
    if rand() < external_total/(external_total + internal_total)
      external[i] = true
    else
      source = findfirst(rand(Multinomial(1, internal_rates[:, i]/internal_total)))
      internal[source, i] = true
    end
  end
  return Network(external, internal)
end


"""
propose(riskpriors::SEIR_RiskParameterPriors)

Randomly generate a set of `SEIR_RiskParameters` from their prior distributions
"""
function propose(riskpriors::SEIR_RiskParameterPriors)
  sparks = Float64[]
  susceptibility = Float64[]
  transmissibility = Float64[]
  infectivity = Float64[]
  latency = Float64[]
  removal = Float64[]

  for i = 1:length(riskpriors.sparks)
    push!(sparks, rand(riskpriors.sparks[i]))
  end

  for i = 1:length(riskpriors.susceptibility)
    push!(susceptibility, rand(riskpriors.susceptibility[i]))
  end

  for i = 1:length(riskpriors.transmissibility)
    push!(transmissibility, rand(riskpriors.transmissibility[i]))
  end

  for i = 1:length(riskpriors.infectivity)
    push!(infectivity, rand(riskpriors.infectivity[i]))
  end

  for i = 1:length(riskpriors.latency)
    push!(latency, rand(riskpriors.latency[i]))
  end

  for i = 1:length(riskpriors.removal)
    push!(removal, rand(riskpriors.removal[i]))
  end

  return SEIR_RiskParameters(sparks,
                             susceptibility,
                             transmissibility,
                             infectivity,
                             latency,
                             removal)
end


"""
propose(riskpriors::SIR_RiskParameterPriors)

Randomly generate a set of `SIR_RiskParameters` from their prior distributions
"""
function propose(riskpriors::SIR_RiskParameterPriors)
  sparks = Float64[]
  susceptibility = Float64[]
  transmissibility = Float64[]
  infectivity = Float64[]
  removal = Float64[]

  for i = 1:length(riskpriors.sparks)
    push!(sparks, rand(riskpriors.sparks[i]))
  end

  for i = 1:length(riskpriors.susceptibility)
    push!(susceptibility, rand(riskpriors.susceptibility[i]))
  end

  for i = 1:length(riskpriors.transmissibility)
    push!(transmissibility, rand(riskpriors.transmissibility[i]))
  end

  for i = 1:length(riskpriors.infectivity)
    push!(infectivity, rand(riskpriors.infectivity[i]))
  end

  for i = 1:length(riskpriors.removal)
    push!(removal, rand(riskpriors.removal[i]))
  end

  return SIR_RiskParameters(sparks,
                            susceptibility,
                            transmissibility,
                            infectivity,
                            removal)
end


"""
propose(riskpriors::SEI_RiskParameterPriors)

Randomly generate a set of `SEI_RiskParameters` from their prior distributions
"""
function propose(riskpriors::SEI_RiskParameterPriors)
  sparks = Float64[]
  susceptibility = Float64[]
  transmissibility = Float64[]
  infectivity = Float64[]
  latency = Float64[]

  for i = 1:length(riskpriors.sparks)
    push!(sparks, rand(riskpriors.sparks[i]))
  end

  for i = 1:length(riskpriors.susceptibility)
    push!(susceptibility, rand(riskpriors.susceptibility[i]))
  end

  for i = 1:length(riskpriors.transmissibility)
    push!(transmissibility, rand(riskpriors.transmissibility[i]))
  end

  for i = 1:length(riskpriors.infectivity)
    push!(infectivity, rand(riskpriors.infectivity[i]))
  end

  for i = 1:length(riskpriors.latency)
    push!(latency, rand(riskpriors.latency[i]))
  end

  return SEI_RiskParameters(sparks,
                            susceptibility,
                            transmissibility,
                            infectivity,
                            latency)
end


"""
propose(riskpriors::SI_RiskParameterPriors)

Randomly generate a set of `SI_RiskParameters` from their prior distributions
"""
function propose(riskpriors::SIR_RiskParameterPriors)
  sparks = Float64[]
  susceptibility = Float64[]
  transmissibility = Float64[]
  infectivity = Float64[]

  for i = 1:length(riskpriors.sparks)
    push!(sparks, rand(riskpriors.sparks[i]))
  end

  for i = 1:length(riskpriors.susceptibility)
    push!(susceptibility, rand(riskpriors.susceptibility[i]))
  end

  for i = 1:length(riskpriors.transmissibility)
    push!(transmissibility, rand(riskpriors.transmissibility[i]))
  end

  for i = 1:length(riskpriors.infectivity)
    push!(infectivity, rand(riskpriors.infectivity[i]))
  end

  return SI_RiskParameters(sparks,
                           susceptibility,
                           transmissibility,
                           infectivity)
end


"""
propose(currentstate::SEIR_RiskParameters,
        variance::Array{Float64, 2})

Generate a `SEIR_RiskParameters` proposal using the multivariate normal distribution
as the transition kernel, with a previous set of `SEIR_RiskParameters` as the mean
vector
"""
function propose(currentstate::SEIR_RiskParameters,
                 variance::Array{Float64, 2})
  newstate = rand(MvNormal(Vector(currentstate), variance))
  inds = cumsum([length(currentstate.sparks);
                 length(currentstate.susceptibility);
                 length(currentstate.transmissibility);
                 length(currentstate.infectivity);
                 length(currentstate.latency);
                 length(currentstate.removal)])
  return SEIR_RiskParameters(newstate[1:inds[1]],
                             newstate[inds[1]+1:inds[2]],
                             newstate[inds[2]+1:inds[3]],
                             newstate[inds[3]+1:inds[4]],
                             newstate[inds[4]+1:inds[5]],
                             newstate[inds[5]+1:inds[6]])
end


"""
propose(currentstate::SIR_RiskParameters,
        variance::Array{Float64, 2})

Generate a `SIR_RiskParameters` proposal using the multivariate normal distribution
as the transition kernel, with a previous set of `SIR_RiskParameters` as the mean
vector
"""
function propose(currentstate::SIR_RiskParameters,
                 variance::Array{Float64, 2})
  newstate = rand(MvNormal(Vector(currentstate), variance))
  inds = cumsum([length(currentstate.sparks);
                 length(currentstate.susceptibility);
                 length(currentstate.transmissibility);
                 length(currentstate.infectivity);
                 length(currentstate.removal)])
  return SIR_RiskParameters(newstate[1:inds[1]],
                            newstate[inds[1]+1:inds[2]],
                            newstate[inds[2]+1:inds[3]],
                            newstate[inds[3]+1:inds[4]],
                            newstate[inds[4]+1:inds[5]])
end


"""
propose(currentstate::SEI_RiskParameters,
        variance::Array{Float64, 2})

Generate a `SEI_RiskParameters` proposal using the multivariate normal distribution
as the transition kernel, with a previous set of `SEI_RiskParameters` as the mean
vector
"""
function propose(currentstate::SEI_RiskParameters,
                 variance::Array{Float64, 2})
  newstate = rand(MvNormal(Vector(currentstate), variance))
  inds = cumsum([length(currentstate.sparks);
                 length(currentstate.susceptibility);
                 length(currentstate.transmissibility);
                 length(currentstate.infectivity);
                 length(currentstate.latency)])
  return SEI_RiskParameters(newstate[1:inds[1]],
                            newstate[inds[1]+1:inds[2]],
                            newstate[inds[2]+1:inds[3]],
                            newstate[inds[3]+1:inds[4]],
                            newstate[inds[4]+1:inds[5]])
end


"""
propose(currentstate::SI_RiskParameters,
        variance::Array{Float64, 2})

Generate a `SI_RiskParameters` proposal using the multivariate normal distribution
as the transition kernel, with a previous set of `SI_RiskParameters` as the mean
vector
"""
function propose(currentstate::SI_RiskParameters,
                 variance::Array{Float64, 2})
  newstate = rand(MvNormal(Vector(currentstate), variance))
  inds = cumsum([length(currentstate.sparks);
                 length(currentstate.susceptibility);
                 length(currentstate.transmissibility);
                 length(currentstate.infectivity)])
  return SI_RiskParameters(newstate[1:inds[1]],
                           newstate[inds[1]+1:inds[2]],
                           newstate[inds[2]+1:inds[3]],
                           newstate[inds[3]+1:inds[4]])
end


"""
propose(i::Int64,
        j::Int64,
        events::SEIR_Events,
        network::Network,
        observations::EventObservations,
        variance::Float64,
        exposureextent::Float64,
        infectionextent::Float64,
        removalextent::Float64)

Event time augmentation for a single event time
"""
function propose(i::Int64,
                 j::Int64,
                 events::SEIR_Events,
                 network::Network,
                 observations::EventObservations,
                 variance::Float64,
                 exposureextent::Float64,
                 infectionextent::Float64,
                 removalextent::Float64)
  exposed = copy(events.exposed)
  infected = copy(events.infected)
  removed = copy(events.removed)
  pathfrom = pathwayfrom(i, network, 1)
  pathto = pathwayto(i, network, 2)
  # Exposure time
  if j == 1
    if length(pathto) > 1
      exposure_lb = maximum([infected[i]-exposureextent;
                             infected[pathto[2]]])
      exposure_ub = minimum([infected[i];
                             removed[pathto[2]]])
      exposed[i] = rand(TruncatedNormal(exposed[i], variance, exposure_lb, exposure_ub))
    else
      exposure_lb = infected[i]-exposureextent
      exposure_ub = infected[i]
      exposed[i] = rand(TruncatedNormal(exposed[i], variance, exposure_lb, exposure_ub))
    end
  # Infection time
  elseif j == 2
    infection_lb = maximum([exposed[i];
                            observations.infected[i]-infectionextent])
    infection_ub = minimum([exposed[pathfrom[2:end]];
                            observations.infected[i]])
    infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
  # Removal time
  elseif j == 3
    removal_lb = maximum([exposed[pathfrom[2:end]];
                          observations.infected[i];
                          observations.removed[i]-removalextent])
    removal_ub = observations.removed[i]
    removed[i] = rand(TruncatedNormal(removed[i], variance, removal_lb, removal_ub))
  end
  return SEIR_Events(exposed,
                     infected,
                     removed)
end


"""
propose(i::Int64,
        j::Int64,
        events::SIR_Events,
        network::Network,
        observations::EventObservations,
        variance::Float64,
        infectionextent::Float64,
        removalextent::Float64)

Event time augmentation for a single event time
"""
function propose(i::Int64,
                 j::Int64,
                 events::SIR_Events,
                 network::Network,
                 observations::EventObservations,
                 variance::Float64,
                 infectionextent::Float64,
                 removalextent::Float64)
  infected = copy(events.infected)
  removed = copy(events.removed)
  pathfrom = pathwayfrom(i, network, 1)
  pathto = pathwayto(i, network, 2)
  # Infection time
  elseif j == 1
    if length(pathto) > 1
      infection_lb = maximum([infected[pathto[2]]];
                              observations.infected[i]-infectionextent])
      infection_ub = minimum([removed[pathto[2]]];
                              infected[pathfrom[2:end]];
                              observations.infected[i]])
      infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
    else
      infection_lb = observations.infected[i]-infectionextent
      infection_ub = minimum([infected[pathfrom[2:end]];
                              observations.infected[i]])
      infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
    end
  # Removal time
  elseif j == 2
    removal_lb = maximum([infected[pathfrom[2:end]];
                          observations.infected[i];
                          observations.removed[i]-removalextent])
    removal_ub = observations.removed[i]
    removed[i] = rand(TruncatedNormal(removed[i], variance, removal_lb, removal_ub))
  end
  return SIR_Events(infected,
                    removed)
end


"""
propose(i::Int64,
        j::Int64,
        events::SEI_Events,
        network::Network,
        observations::EventObservations,
        variance::Float64,
        exposureextent::Float64,
        infectionextent::Float64,
        removalextent::Float64)

Event time augmentation for a single event time
"""
function propose(i::Int64,
                 j::Int64,
                 events::SEI_Events,
                 network::Network,
                 observations::EventObservations,
                 variance::Float64,
                 exposureextent::Float64,
                 infectionextent::Float64)
  exposed = copy(events.exposed)
  infected = copy(events.infected)
  pathfrom = pathwayfrom(i, network, 1)
  pathto = pathwayto(i, network, 2)
  # Exposure time
  if j == 1
    if length(pathto) > 1
      exposure_lb = maximum([infected[i]-exposureextent;
                             infected[pathto[2]]])
      exposure_ub = infected[i]
      exposed[i] = rand(TruncatedNormal(exposed[i], variance, exposure_lb, exposure_ub))
    else
      exposure_lb = infected[i]-exposureextent
      exposure_ub = infected[i]
      exposed[i] = rand(TruncatedNormal(exposed[i], variance, exposure_lb, exposure_ub))
    end
  # Infection time
  elseif j == 2
    infection_lb = maximum([exposed[i];
                            observations.infected[i]-infectionextent])
    infection_ub = minimum([exposed[pathfrom[2:end]];
                            observations.infected[i]])
    infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
  end
  return SEI_Events(exposed,
                    infected)
end


"""
propose(i::Int64,
        j::Int64,
        events::SI_Events,
        network::Network,
        observations::EventObservations,
        variance::Float64,
        infectionextent::Float64,
        removalextent::Float64)

Event time augmentation for a single event time
"""
function propose(i::Int64,
                 j::Int64,
                 events::SI_Events,
                 network::Network,
                 observations::EventObservations,
                 variance::Float64,
                 infectionextent::Float64)
  infected = copy(events.infected)
  removed = copy(events.removed)
  pathfrom = pathwayfrom(i, network, 1)
  pathto = pathwayto(i, network, 2)
  # Infection time
  elseif j == 1
    if length(pathto) > 1
      infection_lb = maximum([infected[pathto[2]]];
                              observations.infected[i]-infectionextent])
      infection_ub = minimum([infected[pathfrom[2:end]];
                              observations.infected[i]])
      infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
    else
      infection_lb = observations.infected[i]-infectionextent
      infection_ub = minimum([infected[pathfrom[2:end]];
                              observations.infected[i]])
      infected[i] = rand(TruncatedNormal(infected[i], variance, infection_lb, infection_ub))
    end
  end
  return SI_Events(infected)
end
