struct EventExtents{T <: EpidemicModel}
  exposure::Float64
  infection::Float64
  removal::Float64

  EventExtents{SEIR}(exposure, infection, removal) = new(exposure, infection, removal)
  EventExtents{SEI}(exposure, infection) = new(exposure, infection, NaN)
  EventExtents{SIR}(infection, removal) = new(NaN, infection, removal)
  EventExtents{SI}(infection) = new(NaN, infection, NaN)
end
