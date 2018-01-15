abstract type Event end

type SEIR_Event <: Event
  eventtype::Int64
  eventindex::Int64

  function SEIR_Event(eventtype::Int64,
                      eventindex::Int64)
    if !(0 <= eventtype <= 4)
      error("Invalid event type")
    end
    return new(eventtype, eventindex)
  end
end

type SIR_Event <: Event
  eventtype::Int64
  eventindex::Int64

  function SIR_Event(eventtype::Int64,
                     eventindex::Int64)
    if !(0 <= eventtype <= 3)
      error("Invalid event type")
    end
    return new(eventtype, eventindex)
  end
end

type SEI_Event <: Event
  eventtype::Int64
  eventindex::Int64

  function SEI_Event(eventtype::Int64,
                     eventindex::Int64)
    if !(0 <= eventtype <= 3)
      error("Invalid event type")
    end
    return new(eventtype, eventindex)
  end
end

type SI_Event <: Event
  eventtype::Int64
  eventindex::Int64

  function SI_Event(eventtype::Int64,
                    eventindex::Int64)
    if !(0 <= eventtype <= 2)
      error("Invalid event type")
    end
    return new(eventtype, eventindex)
  end
end


function getindex(x::Event, i::Int64)
  if i == 1
    x.eventtype
  elseif i == 2
    x.eventindex
  else
    throw(BoundsError)
  end
end
