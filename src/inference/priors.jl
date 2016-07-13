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
