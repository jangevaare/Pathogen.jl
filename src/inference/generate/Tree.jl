function generate(::Type{Tree},
                  events::Events{S},
                  obs::EventObservations{S, M},
                  network::TransmissionNetwork;
                  mrca::Float64=0.0) where {
                  S <: DiseaseStateSequence,
                  M <: ILM}
  return generate(Tree, events, obs.infection, network, mrca=mrca)
end