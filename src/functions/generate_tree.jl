"""
generate_tree(events::Events,
              observations::EventObservations,
              network::Network;
              roottime=0.::Float64)

Generates a phylogenetic tree based on observations and transmission network
"""
function generate_tree(events::Events,
                       observations::EventObservations,
                       network::Network;
                       roottime=0.::Float64)
  # Initialization and error checking
  # Determine significant exposures (or infections)
  # Exposures that result in a observation are considered significant
  events_type = typeof(events)
  tree = Tree()
  pathways = [pathwayfrom(i, network) for i = 1:events.individuals]
  eventtimes = fill(NaN, (events.individuals, 2))
  significant_exposures = Int64[]
  if events_type == SEIR_Events
    if typeof(observations) != SEIR_EventObservations
      error("Mismatch in the subtypes of `Events` and `EventObservations`")
    end
    eventtimes = [events.exposed observations.infected]
    significant_exposures = [any(observations.infected[pathways[i]] .>= events.exposed[i]) for i = 1:events.individuals]
  elseif events_type == SIR_Events
    if typeof(observations) != SIR_EventObservations
      error("Mismatch in the subtypes of `Events` and `EventObservations`")
    end
    eventtimes = [events.infected observations.infected]
    significant_exposures = [any(observations.infected[pathways[i]] .>= events.infected[i]) for i = 1:events.individuals]
  elseif events_type == SEI_Events
    if typeof(observations) != SEI_EventObservations
      error("Mismatch in the subtypes of `Events` and `EventObservations`")
    end
    eventtimes = [events.exposed observations.infected]
    significant_exposures = [any(observations.infected[pathways[i]] .>= events.exposed[i]) for i = 1:events.individuals]
  elseif events_type == SI_Events
    if typeof(observations) != SI_EventObservations
      error("Mismatch in the subtypes of `Events` and `EventObservations`")
    end
    eventtimes = [events.infected observations.infected]
    significant_exposures = [any(observations.infected[pathways[i]] .>= events.infected[i]) for i = 1:events.individuals]
  end
  eventorder = sortperm(eventtimes[:])
  eventnodes = fill(Nullable{Int64}(), (events.individuals, 3))

  # Add a node for each infection observation event
  infection_observations = find(.!isnan.(observations.infected))
  addnodes!(tree, length(infection_observations))

  for i = 1:length(infection_observations)
    eventnodes[infection_observations[i], 3] = i
  end

  # Add a root node
  addnode!(tree)
  rootnode = length(tree.nodes)

  # Iterate through all events to build tree
  for i = 1:length(eventorder)

    # Stop when no remaining valid event times
    isnan(eventtimes[eventorder[i]]) && break

    # Determine the individual and event type
    event = ind2sub(size(eventtimes), eventorder[i])

    # For exposure events...
    if event[2] == 1

      # Add exposure event to tree only if it is significant
      if significant_exposures[event[1]]

        # For significant external exposure events...
        if network.external[event[1]]

          # Parent node of external exposure is the root of the tree...
          parentnode = rootnode
          branch_length = eventtimes[eventorder[i]] - roottime

        # For significant internal exposure events...
        else

          # Determine exposure source
          source = findfirst(network.internal[:, event[1]])

          # Previous significant exposures from this exposure source
          source_priorexposures = network.internal[source, :][:] .& (eventtimes[:, 1] .< eventtimes[eventorder[i]]) .& significant_exposures

          # For when the exposure source has yet to be detected...
          if isnan(eventtimes[source, 2]) || eventtimes[source, 2] > eventtimes[event[1], 1]

            # When there have been prior exposures from this undetected source
            if any(source_priorexposures)
              parentnode = get(eventnodes[source_priorexposures, 1][indmax(eventtimes[source_priorexposures, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorexposures, 1])

            # When there have not been any prior exposures from this undetected source
            else
              parentnode = get(eventnodes[source, 1])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 1]
            end

          # For when the exposure source has been detected...
          else
            # And detection of exposure source is most recent, relevant event...
            if !any(source_priorexposures) || all(eventtimes[source, 2] .> eventtimes[source_priorexposures, 1])
              parentnode = get(eventnodes[source, 2])
              branch_length = eventtimes[eventorder[i]] - eventtimes[source, 2]

            # Or some prior exposure is most recent, relevant event...
            else
              parentnode = get(eventnodes[source_priorexposures, 1][indmax(eventtimes[source_priorexposures, 1])])
              branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[source_priorexposures, 1])
            end
          end
        end

        # Add branch and node to tree
        branch!(tree, parentnode, branch_length)
        eventnodes[eventorder[i]] = length(tree.nodes)
      end

    # Infection observation event
    elseif event[2] == 2

      # Individual has had significant exposures prior to being observed as infected
      priorexposures = network.internal[event[1], :][:] .& (eventtimes[:, 1] .< eventtimes[eventorder[i]]) .& significant_exposures
      if any(priorexposures)
        parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
        branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])

      # Individual has not exposed others prior to being observed as infected
      else
        parentnode = get(eventnodes[event[1], 1])
        branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
      end

      # Individual goes on to have significant exposures
      if any(network.internal[event[1], :][:] .& (eventtimes[:, 1] .>= eventtimes[eventorder[i]]) .& significant_exposures)
        branch!(tree, parentnode, branch_length)
        eventnodes[event[1], 2] = length(tree.nodes)

        # Add zero-length branch so observation event is a leaf node
        addbranch!(tree, get(eventnodes[event[1], 2]), get(eventnodes[event[1], 3]), 0.)

      # Individual does not go one to have any significant exposures
      else
        addbranch!(tree, parentnode, get(eventnodes[event[1], 3]), branch_length)
      end
    end
  end
  return tree
end
