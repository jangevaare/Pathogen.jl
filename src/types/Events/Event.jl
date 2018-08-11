mutable struct Event{T <: EpidemicModel}
  time::Float64
  individual::Int64
  new_state::DiseaseState

  function Event{T}(time::Float64) where T <: EpidemicModel
    if time !== Inf
      @error "This initialization method is intended for when for generated event times of Inf, i.e. when an epidemic is complete and no further events are possible."
    end
    x = new{T}()
    x.time = Inf
    return x
  end

  function Event{T}(time::Float64, id::Int64, new_state::DiseaseState) where T <: EpidemicModel
    if new_state in _state_progressions[T]
      return new{T}(time, id, new_state)
    else
      @error "Invalid disease state provided for $(T) Epidemic Models."
    end
  end

  function Event(new_time::Float64, event::Event{T}) where T <: EpidemicModel
    return Event{T}(new_time, event.individual, event.new_state)
  end
end

function Base.copy(x::Event{T}) where T <: EpidemicModel
  return Event{T}(copy(x.time), copy(x.id), copy(x.new_state))
end

function _new_transmission(e::Event{T}) where T <: Union{SEIR, SEI}
  return e.new_state == State_E
end

function _new_transmission(e::Event{T}) where T <: Union{SIR, SI}
  return e.new_state == State_I
end

function _time(x::Event{T}) where T <: EpidemicModel
  return x.time
end

function _individual(x::Event{T}) where T <: EpidemicModel
  return x.individual
end

function _new_state(x::Event{T}) where T <: EpidemicModel
  return x.new_state
end

function Base.show(io::IO, x::Event{T}) where T <: EpidemicModel
  return print(io, "Transition of individual $(x.individual) into $(x.new_state) state at t = $(x.time)")
end
