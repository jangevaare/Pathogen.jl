abstract EventExtents

type SEIR_EventExtents <: EventExtents
  exposure::Float64
  infection::Float64
  removal::Float64
end


type SIR_EventExtents <: EventExtents
  infection::Float64
  removal::Float64
end


type SEI_EventExtents <: EventExtents
  exposure::Float64
  infection::Float64
end


type SI_EventExtents <: EventExtents
  infection::Float64
end
