function _accept(rng::AbstractRNG, lp1::Float64,
                 lp2::Float64)
  return rand(rng) <= exp(lp1 - lp2)
end


_accept(lp1::Float64, lp2::Float64) = _accept(_default_rng(), lp1, lp2)
