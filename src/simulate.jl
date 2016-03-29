"""
Stores rate information for simulation purposes
"""
type Rates
  internal_exposure::Array{Float64,2}
  external_exposure::Vector{Float64}
  infection::Vector{Float64}
  detection::Vector{Float64}
  removal::Vector{Float64}

  function Rates(individuals::Int64)
    if individuals <= 0
      error("At least one individual must be present")
    end
    return new(fill(0., (individuals, individuals)),
               fill(0., individuals),
               fill(0., individuals),
               fill(0., individuals),
               fill(0., individuals))
  end
end


"""
Simulate events and a transmission tree for a
phylodynamic individual level model of infectious disease
"""
function simulate(n::Int64,
                  population::DataFrame,
                  risk_funcs::RiskFunctions,
                  risk_params::RiskParameters,
                  index_case=1::Int64)

  individuals = size(DataFrame, 2)

  # Initialize rate arrays
  rates = Rates(individuals)
  rates.infection[index_case] = risk_funcs.latency([risk_params.latency], population, index_case)
  for i in [1:(index_case-1); (index_case+1):individuals]
    rates.external_exposure[i] = risk_funcs.sparks([risk_params.sparks], population, i)
  end

  # Initialize events data frame
  events = DataFrame(time=Float64[], event=Tuple[], node=Tuple[])
  push!(events, [0. (2, index_case) (1, 1)])

  # Initialize tree vector
  trees = Tree[]
  push!(trees, Tree())
  add_node!(trees[1])

  info("Initialization complete")



  # rates.internal_exposure[k, i] = risk_funcs.susceptibility([risk_params.susceptibility], population, i) *
  #                                          risk_funcs.transmissibility([risk_params.transmissibility], population, k) *
  #                                          risk_funcs.infectivity([risk_params.infectivity], population, i, k)

end
