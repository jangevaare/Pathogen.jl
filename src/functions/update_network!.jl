"""
update_network!(network::Network,
                event::Event

A function to update a `Network` object based on an event occurence
"""
function update_network!(network::Network,
                         event::Event)
  if event[1] == 1
    network.external[event[2]] = true
  elseif event[1] == 2
    network.internal[event[2]] = true
  end
  return network
end
