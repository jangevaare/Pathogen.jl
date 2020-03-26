function generate(::Type{PhyloTree},
                  events::Events{S},
                  obs::Vector{Float64},
                  network::TransmissionNetwork,
                  mrca::Float64) where {
                  S <: DiseaseStateSequence}

  if length(obs) != events.individuals
    throw(ErrorException("Infection observation vector does not match number of individuals in population"))
  elseif mrca > events.start_time
    throw(ErrorException("Most recent common ancestor to external transmissions must be <= the start of an epidemic"))
  end

  # Initialization
  # trees = Tree[]
  tree = Tree(mrca)
  addnode!(tree) # Root/MRCA node
  local event_times

  if S in [SEIR; SEI]
    event_times = [events.exposure obs]
  elseif S in [SIR; SI]
    event_times = [events.infection obs]
  else
    throw(ErrorException("Unrecognized DiseaseStateSequence"))
  end

  event_nodes = Array{Union{Nothing, Int64}, 2}(fill(nothing, size(event_times)))
  event_order = sortperm(event_times[:])
  event_lookup = CartesianIndices(size(event_times))

  pathways = [_pathway_from(i, network) for i = 1:events.individuals]

  # Determine significant transmissions (results in a observation in the future)
  significant_transmissions = [event_times[i, 1] == -Inf || any(event_times[pathways[i], 2] .>= Ref(event_times[i, 1])) for i = 1:events.individuals]

  # Iterate through all events to build tree
  for i = eachindex(event_order)

    # Stop when no remaining valid event times
    isnan(event_times[event_order[i]]) && break

    # Determine the individual and event type
    event_individual, event_type = Tuple(event_lookup[event_order[i]])

    @debug "Building phylogenetic tree..." Event = i Individual = event_individual Type = event_type == 1 ? "Transmission" : "Observation"

    # For transmission events...
    if event_type == 1

      # Add transmission event to tree only if it is significant
      if significant_transmissions[event_individual]

        # For significant external transmission events...
        if network.external[event_individual]
          @debug "Event $i is a significant external transmission"
          parent_node = 1
          branch_length = mrca - event_times[event_order[i]]
          branch!(tree, parent_node, branch_length)

          # Record tree and node ID
          event_nodes[event_order[i]] = length(tree.nodes)

        # For other significant transmission events...
        else

          # Determine source of transmission
          source = findfirst(network.internal[:, event_individual])

          # If infected initial condition, no source
          if source == nothing
            @debug "Event $i is an initial condition"

            parent_node = 1
            branch_length = mrca - events.start_time

            # Add branch and node
            branch!(tree, parent_node, branch_length)

            # Record node id
            event_nodes[event_order[i]] = length(tree.nodes)

          # Previous significant transmissions from this source of transmission
          else
            @debug "Event $i is a significant internal transmission"

            source_prior_transmissions = findall(network.internal[source, :][:] .& (event_times[:, 1] .< event_times[event_order[i]]) .& significant_transmissions)

            # For when the source of transmission has yet to be detected...
            if isnan(event_times[source, 2]) || event_times[source, 2] > event_times[event_individual, 1]

              # When there have been prior transmissions from this undetected source
              if length(source_prior_transmissions) > 0
                parent_node = event_nodes[source_prior_transmissions[argmax(event_times[source_prior_transmissions, 1])], 1]
                branch_length = event_times[event_order[i]] - maximum(event_times[source_prior_transmissions, 1])

              # When there have not been any prior transmissions from this undetected source
              else
                parent_node = event_nodes[source, 1]
                branch_length = event_times[event_order[i]] - event_times[source, 1]
              end

            # For when the source of transmission has been detected...
            else
              # And detection of source of transmission is most recent, relevant event...
              if length(source_prior_transmissions) == 0 || all(event_times[source, 2] .> event_times[source_prior_transmissions, 1])

                # Determine parent node of the observation leaf node
                parent_node = parentnode(tree, event_nodes[source, 2])
                branch_length = event_times[event_order[i]] - event_times[source, 2]

              # Or some prior transmission is most recent, relevant event...
              else
                parent_node = event_nodes[source_prior_transmissions[argmax(event_times[source_prior_transmissions, 1])], 1]
                branch_length = event_times[event_order[i]] - maximum(event_times[source_prior_transmissions, 1])
              end
            end
            # Add branch and node
            branch!(tree, parent_node, branch_length)

            # Record node id
            event_nodes[event_order[i]] = length(tree.nodes)
          end
        end
      else
        @debug "Event $i is an insignificant transmission"
      end

    # Infection observation event
    elseif event_type == 2

      # Determine if individual has had significant transmissions prior to being observed as infected
      prior_transmissions = findall(network.internal[event_individual, :][:] .&
                            (event_times[:, 1] .< event_times[event_order[i]]) .&
                            significant_transmissions)

      # Parent node is the internal node associated with the final pre-observation transmission
      if length(prior_transmissions) > 0
        parent_node = event_nodes[prior_transmissions[argmax(event_times[prior_transmissions, 1])], 1]
        branch_length = event_times[event_order[i]] - maximum(event_times[prior_transmissions, 1])

      # Individual has not transmitted to others prior to being observed as infected
      else
        parent_node = event_nodes[event_individual, 1]
        branch_length = event_times[event_order[i]] - event_times[event_individual, 1]
      end

      # Individual goes on to have significant transmissions
      if any(network.internal[event_individual, :][:] .& (event_times[:, 1] .>= event_times[event_order[i]]) .& significant_transmissions)

        # Internal node for observation
        branch!(tree, parent_node, branch_length)

        # Add zero-length branch so observation event will be a leaf node
        branch!(tree, length(tree.nodes), 0.)

      # Individual does not go one to have any significant transmissions
      # so will surely be a leaf node
      else
        branch!(tree, parent_node, branch_length)
      end
      # Record leaf node id
      event_nodes[event_order[i]] = length(tree.nodes)
    end
  end
  return tree, event_nodes[:,2]
end