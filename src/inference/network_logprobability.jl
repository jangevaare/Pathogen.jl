# These are the same source weights used by generate(Transmission, ...).
function _transmission_weight(tr::TransmissionRates, prior, source::Integer, id::Integer)
  rate = source == 0 ? tr.external[id] : tr.internal[source, id]
  prior === nothing && return rate
  return rate * (source == 0 ? prior.external[id] : prior.internal[source, id])
end

function _transmission_total(tr::TransmissionRates, prior, id::Integer)
  total = 0.0
  for source in 0:individuals(tr)
    weight = _transmission_weight(tr, prior, source, id)
    isfinite(weight) && weight >= 0.0 || return NaN
    total += weight
  end
  return total
end

function _has_network_support(tr::TransmissionRates, prior, transmissions)
  # A finite epidemic likelihood already establishes support without a prior.
  prior === nothing && return true
  return all(id -> 0.0 < _transmission_total(tr, prior, id) < Inf, transmissions)
end

function _network_logprobabilities(network::TransmissionNetwork, tr::TransmissionRates,
                                   prior, transmissions, starting_states)
  n = individuals(tr)
  length(starting_states) == individuals(network) == n ||
    throw(DimensionMismatch("Network, rates and initial states must have equal lengths"))
  terms = zeros(n)
  occurred = falses(n)
  occurred[transmissions] .= true
  for id in 1:n
    source = 0
    incoming = Int(network.external[id])
    for k in nzrange(network.internal, id)
      if nonzeros(network.internal)[k]
        incoming += 1
        source = rowvals(network.internal)[k]
      end
    end
    if !occurred[id]
      expected_external = starting_states[id] != State_S
      if incoming != Int(expected_external) || network.external[id] != expected_external
        terms[id] = -Inf
      end
      continue
    end
    if incoming != 1 || starting_states[id] != State_S
      terms[id] = -Inf
      continue
    end
    selected = _transmission_weight(tr, prior, source, id)
    total = _transmission_total(tr, prior, id)
    terms[id] = selected > 0.0 && isfinite(total) && total > 0.0 ?
                log(selected) - log(total) : -Inf
  end
  return terms
end

# Redrawn parents are proposed from their source distribution. Their target
# probabilities cancel against that proposal ratio; other recipients still count.
function _retained_network_logratio(proposed, current, redrawn)
  skipped = falses(length(current))
  skipped[redrawn] .= true
  delta = 0.0
  for id in eachindex(current)
    skipped[id] || (delta += proposed[id] - current[id])
  end
  return delta
end
