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
  events = DataFrame(time=Float64[], event=Tuple[], node=Tuple[])
  trees = Tree[]
  time = 0.0

  # Initialize rate arrays
  rates = Rates(size(DataFrame, 2))
end
