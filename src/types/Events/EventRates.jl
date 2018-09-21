mutable struct EventRates{T <: EpidemicModel}
  exposure::Vector{Float64}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64

  function EventRates{T}(n::Int64) where T <: SEIR
    x = new{T}(fill(0.0, n), fill(0.0, n), fill(0.0, n), n)
    return x
  end

  function EventRates{T}(n::Int64) where T <: SEI
    x = new{T}()
    x.exposure = fill(0.0, n)
    x.infection = fill(0.0, n)
    x.individuals = n
    return x
  end

  function EventRates{T}(n::Int64) where T <: SIR
    x = new{T}()
    x.infection = fill(0.0, n)
    x.removal = fill(0.0, n)
    x.individuals = n
    return x
  end

  function EventRates{T}(n::Int64) where T <: SI
    x = new{T}()
    x.infection = fill(0.0, n)
    x.individuals = n
    return x
  end
end

function Base.getindex(x::EventRates{T}, new_state::DiseaseState) where T <: EpidemicModel
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    @error "Unrecognized indexing disease state"
  end
end

function Base.getindex(x::EventRates{T}, new_states::Vector{DiseaseState}) where T <: EpidemicModel
  y = x[new_states[1]]
  for i=2:length(new_states)
    y = hcat(y, x[new_states[i]])
  end
  return y
end
