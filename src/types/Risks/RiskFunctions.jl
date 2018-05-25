mutable struct RiskFunctions{T <: EpidemicModel}
  sparks::Function
  susceptibility::Function
  transmissibility::Function
  infectivity::Function
  latency::Function
  removal::Function

  function RiskFunctions{T}(sp, su, tr, in, la, re) where {T <: SEIR}

    return new{T}(sp, su, tr, in, la, re)
  end

  function RiskFunctions{T}(sp, su, tr, in, la) where {T <: SEI}

    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    x.latency = la
    return x
  end

  function RiskFunctions{T}(sp, su, tr, in, re) where {T <: SIR}

    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    x.removal = re
    return x
  end

  function RiskFunctions{T}(sp, su, tr, in) where {T <: SI}

    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.transmissibility = tr
    x.infectivity = in
    return x
  end
end
