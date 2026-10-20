# Scratch state belongs to a single serial likelihood stream. Returned rate and
# transmission buffers are borrowed until the next evaluation with this workspace.
mutable struct LikelihoodWorkspace{T <: DiseaseStateSequence}
  disease_states::DiseaseStates
  transmission_rates::TransmissionRates
  event_rates::EventRates{T}
  event_times::Vector{Float64}
  event_order::Vector{Int}
  transmissions::Vector{Int64}
  rate_total::Base.RefValue{Float64}
  proposal_cache::Union{Nothing, TransmissionRateCache}
end

function LikelihoodWorkspace(::Type{T}, n::Integer) where T <: DiseaseStateSequence
  return LikelihoodWorkspace{T}(fill(State_S, n), TransmissionRates(n),
                                EventRates{T}(Int64(n)), Float64[], Int[], Int64[], Ref(0.0), nothing)
end

function _event_order!(workspace::LikelihoodWorkspace{T}, events::Events{T}) where T
  n = individuals(events)
  states = convert(DiseaseStates, T)[2:end]
  count = n * length(states)
  changed = length(workspace.event_times) != count
  resize!(workspace.event_times, count)
  for (j, state) in enumerate(states), i in 1:n
    k = i + (j-1)*n
    time = events[state][i]
    # Check values as well as identity: callers may edit an Events vector in place.
    if changed || !isequal(workspace.event_times[k], time)
      workspace.event_times[k] = time
      changed = true
    end
  end
  if changed
    resize!(workspace.event_order, count)
    sortperm!(workspace.event_order, workspace.event_times)
  end
  return workspace.event_order
end


function _proposal_cache!(workspace::LikelihoodWorkspace, n::Integer)
  cache = workspace.proposal_cache
  if cache === nothing
    cache = TransmissionRateCache(n)
    workspace.proposal_cache = cache
  else
    fill!(cache.external, nothing)
    fill!(cache.internal, nothing)
  end
  return cache
end
