struct Event{T <: EpidemicModel}
  time::Float64
  individual::Int64
  new_state::DiseaseState

  function Event{T}(time::Float64) where T <: EpidemicModel
    if time !== Inf
      error("This initialization method is intended for when for generated event times of Inf, i.e. when an epidemic is complete and no further events are possible.")
    end
    x = new{T}()
    x.time = Inf
    return x
  end

  function Event{T}(time::Float64, id::Int64, new_state::DiseaseState) where T <: EpidemicModel
    new_state in _state_progressions[T] || error("Invalid disease state provided for $(T) Epidemic Models.")
    return new{T}(time, id, new_state)
  end
end
