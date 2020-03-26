function update!(states::DiseaseStates,
                 event::Event{S}) where S <: DiseaseStateSequence
  states[event.individual] = event.new_state
  return states
end

