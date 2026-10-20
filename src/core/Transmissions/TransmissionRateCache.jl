mutable struct TransmissionRateCache <: AbstractTN
  external::Vector{Union{Nothing, Float64}}
  internal::Array{Union{Nothing, Float64}, 2}

  function TransmissionRateCache(x::Integer)
    return new(Union{Nothing, Float64}[nothing for i = 1:x],
               Union{Nothing, Float64}[nothing for i = 1:x, j = 1:x])
  end
end

# Keep the caller's cache object, but retain the buffers belonging to the accepted
# parameters. Rebinding a local variable in update! does not update iterate!'s cache.
function _retain_cache!(destination::TransmissionRateCache, accepted::TransmissionRateCache)
  destination.external = accepted.external
  destination.internal = accepted.internal
  return destination
end
_retain_cache!(::Nothing, accepted) = nothing
