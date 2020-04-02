abstract type AbstractEvent end

struct NoEvent <: AbstractEvent end

struct Event{S <: DiseaseStateSequence} <: AbstractEvent
  time::Float64
  individual::Int64
  new_state::DiseaseState

  function Event{S}(time::Float64, id::Int64, new_state::DiseaseState) where S <: DiseaseStateSequence
    if new_state in convert(DiseaseStates, S)
      return new{S}(time, id, new_state)
    else
      @error "Invalid disease state provided for $S model."
    end
  end

  function Event(new_time::Float64, event::Event{S}) where S <: DiseaseStateSequence
    return new{S}(new_time, event.individual, event.new_state)
  end
end

function _time(x::Event{S}) where S <: DiseaseStateSequence
  return x.time
end

function _time(x::NoEvent)
  return Inf
end

# function Base.copy(x::Event{S}) where S <: DiseaseStateSequence
#   return Event{S}(copy(x.time), copy(x.id), copy(x.new_state))
# end

function _new_transmission(e::Event{S}) where S <: Union{SEIR, SEI}
  return e.new_state == State_E
end

function _new_transmission(e::Event{S}) where S <: Union{SIR, SI}
  return e.new_state == State_I
end

# function _individual(x::Event{S}) where S <: DiseaseStateSequence
#   return x.individual
# end

# function _new_state(x::Event{S}) where S <: DiseaseStateSequence
#   return x.new_state
# end

function Base.show(io::IO, x::Event{S}) where S <: DiseaseStateSequence
  return print(io, "Transition of individual $(x.individual) into $(x.new_state) state at t = $(x.time)")
end
