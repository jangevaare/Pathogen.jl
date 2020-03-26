
function generate(::Type{Event},
                  last_event::Event{S},
                  σ::Float64,
                  extents::EventExtents{S},
                  obs::EventObservations,
                  events::Events{S}) where {
                  S <: DiseaseStateSequence}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events)
  time = rand(truncated(Normal(_time(last_event), σ),
                        lowerbound,
                        upperbound))
  return Event(time, last_event)
end

function generate(::Type{Event},
                  last_event::Event{S},
                  σ::Float64,
                  extents::EventExtents{S},
                  obs::EventObservations,
                  events::Events{S},
                  network::TransmissionNetwork) where {
                  S <: DiseaseStateSequence}
  lowerbound, upperbound = _bounds(last_event, extents, obs, events, network)
  time = rand(truncated(Normal(_time(last_event), σ),
                        lowerbound,
                        upperbound))
  return Event(time, last_event)
end
