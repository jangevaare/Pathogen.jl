function generate(::Type{TransmissionNetwork},
                  tr::TransmissionRates,
                  tnd::Union{Nothing, TNDistribution},
                  ids::Vector{Int64}) where M <: ILM
  tn = TransmissionNetwork(tr.individuals)                
  for i in ids
    tx = generate(Transmission, tr, tnd, i)
    update!(tn, tx)
  end
  return tn
end

function generate(::Type{TransmissionNetwork},
                  tr::TransmissionRates,
                  mcmc::MCMC,
                  ids::Vector{Int64}) where M <: ILM
  return generate(TransmissionNetwork, tr, mcmc.transmission_network_prior, ids)
end
