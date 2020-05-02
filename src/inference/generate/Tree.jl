function generate(::Type{Tree},
                  events::Events{S},
                  obs::EventObservations{S, M},
                  network::TransmissionNetwork) where {
                  S <: DiseaseStateSequence,
                  M <: ILM}
  return generate(Tree, events, obs.infection, network)
end