mutable struct RiskFunctions{T <: EpidemicModel}
  sparks::Function
  susceptibility::Function
  infectivity::Function
  transmissibility::Function
  latency::Function
  removal::Function

  function RiskFunctions{T}(sp, su, in, tr, la, re) where {T <: SEIR}
    return new{T}(sp, su, in, tr, la, re)
  end

  function RiskFunctions{T}(sp, su, in, tr, la) where {T <: SEI}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.infectivity = in
    x.transmissibility = tr
    x.latency = la
    return x
  end

  function RiskFunctions{T}(sp, su, in, tr, re) where {T <: SIR}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.infectivity = in
    x.transmissibility = tr
    x.removal = re
    return x
  end

  function RiskFunctions{T}(sp, su, in, tr) where {T <: SI}
    x = new{T}()
    x.sparks = sp
    x.susceptibility = su
    x.infectivity = in
    x.transmissibility = tr
    return x
  end
end

function Base.show(io::IO, x::RiskFunctions{T}) where T <: EpidemicModel
  return print(io, "$T model risk functions")
end
