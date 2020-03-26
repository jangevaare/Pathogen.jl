function update!(net::TransmissionNetwork,
                 tx::ExogenousTransmission)
  net.external[tx.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 tx::EndogenousTransmission)
  net.internal[tx.source, tx.individual] = true
  return net
end

function update!(net::TransmissionNetwork,
                 tx::NoTransmission)
  return net
end