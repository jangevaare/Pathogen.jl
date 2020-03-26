function update!(events::Events{S},
                 event::Event{S}) where S <: DiseaseStateSequence
  events[event.new_state][event.individual] = event.time
  return events
end