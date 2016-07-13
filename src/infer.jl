"""
Prior distributions vectors for the `RiskParameters`
"""
type RiskPriors
  sparks::Vector{UnivariateDistribution}
  susceptibility::Vector{UnivariateDistribution}
  transmissibility::Vector{UnivariateDistribution}
  infectivity::Vector{UnivariateDistribution}
  latency::Vector{UnivariateDistribution}
  detection::Vector{UnivariateDistribution}
  removal::Vector{UnivariateDistribution}
end


"""
MCMC trace
"""
type PathogenTrace
  risk_parameters::Vector{RiskParameters}
  events::Vector{Events}
end


"""
MCMC iteration
"""
type PathogenIteration
  risk_parameters::RiskParameters
  events::Events
end


PathogenProposal = PathogenIteration


function push!(trace::PathogenTrace, iteration::PathogenIteration)
  push!(trace.risk_parameters, iteration.risk_parameters)
  push!(trace.events, iteration.events)
  return trace
end


function append!(trace1::PathogenTrace, trace2::PathogenTrace)
  append!(trace1.risk_parameters, trace2.risk_parameters)
  append!(trace1.events, trace2.events)
end


type RiskParameters
  sparks::Vector{Float64}
  susceptibility::Vector{Float64}
  transmissibility::Vector{Float64}
  infectivity::Vector{Float64}
  latency::Vector{Float64}
  detection::Vector{Float64}
  removal::Vector{Float64}
end


function rand(riskpriors::RiskPriors)
  sparks = Float64[]
  susceptibility = Float64[]
  transmissibility = Float64[]
  infectivity = Float64[]
  latency = Float64[]
  detection = Float64[]
  removal = Float64[]

  for i = 1:length(riskpriors.sparks)
    push!(sparks, rand(riskpriors.sparks[i]))
  end

  for i = 1:length(riskpriors.susceptibility)
    push!(susceptibility, rand(riskpriors.susceptibility[i]))
  end

  for i = 1:length(riskpriors.transmissability)
    push!(transmissability, rand(riskpriors.transmissability[i]))
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

  return RiskParameters(sparks,
                        susceptibility,
                        transmissability,
                        infectivity,
                        latency,
                        detection,
                        removal)
end


"""
Event data augmentation
"""
function eventDA(riskparams::RiskParameters, population::DataFrame, detectiontimes::Vector{Float64})
  detectiontimes = fill(NaN, length(detectiontimes))
  for i in 1:length(detectiontimes)
    if !isnan(detectiontimes[i])
      # TODO
    end
  end
  infectiontimes = detectiontimes
  # TODO
end


type Events
  susceptible::Vector{Float64}
  exposed::Vector{Float64}
  infected::Vector{Float64}
  detected::Vector{Float64}
  removed::Vector{Float64}
  network::Vector{Array{Bool}}
