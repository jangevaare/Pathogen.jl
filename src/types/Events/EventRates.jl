mutable struct EventRates{T <: EpidemicModel}
  exposure::Vector{Float64}
  infection::Vector{Float64}
  removal::Vector{Float64}
  individuals::Int64
  # TODO: add reference to `TransmissionRates`?

  function EventRates{T}(individuals::Int64) where T <: EpidemicModel
    return _init_EventRates!(new{T}(), individuals)
  end

  function EventRates{T}(v...) where T <: EpidemicModel
    if unique(length.(v)) != 1
      error("Mismatch in length of rate vectors")
    end
    return _init_EventRates!(new{T}(), v)
  end
end

function _init_EventRates!(x::EventRates{SEIR}, individuals::Int64)
  return _init_EventRates!(x, fill(0., individuals), fill(0., individuals), fill(0., individuals))
end


function _init_EventRates!(x::EventRates{SEI}, individuals::Int64)
  return _init_EventRates!(x, fill(0., individuals), fill(0., individuals))
end


function _init_EventRates!(x::EventRates{SIR}, individuals::Int64)
  return _init_EventRates!(x, fill(0., individuals), fill(0., individuals))
end


function _init_EventRates!(x::EventRates{SI}, individuals::Int64)
  return _init_EventRates!(x, fill(0., individuals))
end

function _init_EventRates!(x::EventRates{T}, v::Tuple) where T <: EpidemicModel
  # Unforeseen consequence of heavy use of incomplete initialization...
  # TODO: Address more elegantly in future
  if length(v) == 3
    return _init_EventRates!(x, v[1], v[2], v[3])
  elseif length(v) == 2
    return _init_EventRates!(x, v[1], v[2])
  elseif length(v) == 1
    return _init_EventRates!(x, v[1])
  else
    error("How did I get here?")
  end
end

function _init_EventRates!(x::EventRates{SEIR}, v...)
  x.exposure = v[1]
  x.infection = v[2]
  x.removal = v[3]
  x.individuals = length(x.infection)
  return x
end

function _init_EventRates!(x::EventRates{SEI}, v...)
  x.exposure = v[1]
  x.infection = v[2]
  x.individuals = length(x.infection)
  return x
end

function _init_EventRates!(x::EventRates{SIR}, v...)
  x.infection = v[1]
  x.removal = v[2]
  x.individuals = length(x.infection)
  return x
end

function _init_EventRates!(x::EventRates{SI}, v...)
  x.infection = v[1]
  x.individuals = length(x.infection)
  return x
end

function Base.getindex(x::EventRates{T}, new_state::DiseaseState) where T <: EpidemicModel
  if new_state == State_E
    return x.exposure
  elseif new_state == State_I
    return x.infection
  elseif new_state == State_R
    return x.removal
  else
    error()
  end
end
