function _accept(lp1::Float64,
                 lp2::Float64)
  return rand() <= exp(lp1 - lp2)
end
