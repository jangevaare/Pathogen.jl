mutable struct Events{T <: EpidemicModel}
  exposure::Vector{Float64}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64

  function Events{T}(e::V, i::V, r::V) where {T <: SEIR, V <: Vector{Float64}}
    if length(unique((length.([e; i; r])))) != 1
      error("Length of event time vectors must be equal")
    end
    return new{T}(e, i, r, length(i))
  end

  function Events{T}(n_ids::Int64) where T <: SEIR
    return Events{T}(fill(NaN, n_ids), fill(NaN, n_ids), fill(NaN, n_ids))
  end

  function Events{T}(e::V, i::V) where {T <: SEI, V <: Vector{Float64}}
    if length(unique((length.([e; i])))) != 1
      error("Length of event time vectors must be equal")
    end
    x = new{T}()
    x.exposure = e
    x.infection = i
    x.individuals = length(i)
    return x
  end

  function Events{T}(n_ids::Int64) where T <: SEI
    return Events{T}(fill(NaN, n_ids), fill(NaN, n_ids))
  end

  function Events{T}(i::V, r::V) where {T <: SIR, V <: Vector{Float64}}
    if length(unique((length.([i; r])))) != 1
      error("Length of event time vectors must be equal")
    end
    x = new{T}()
    x.infection = i
    x.removal = r
    x.individuals = length(i)
    return x
  end

  function Events{T}(n_ids::Int64) where T <: SIR
    return Events{T}(fill(NaN, n_ids), fill(NaN, n_ids))
  end

  function Events{T}(i::V) where {T <: SI, V <: Vector{Float64}}
    x = new{T}()
    x.infection = i
    x.individuals = length(i)
    return x
  end

  function Events{T}(n_ids::Int64) where T <: SI
    return Events{T}(fill(NaN, n_ids))
  end
end

function Base.getindex(x::Events{T}, new_state::DiseaseState) where T <: EpidemicModel
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    error("Unrecognized indexing disease state")
  end
end
