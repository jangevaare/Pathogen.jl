"""
Generate phylogenetic tree based on transmission events
"""
function generatetree(events::Events)
  eventtimes = [events.exposed events.detected events.removed]
  eventorder = sortperm(eventtimes[:])
  eventnodes = fill(Nullable{Int64}(), size(eventtimes))
  tree = Tree()
  for i = 1:length(eventorder)
    isnan(eventtimes[eventorder[i]]) && break
    event = ind2sub(size(eventtimes), eventorder[i])
    if event[2] == 1
      # Exposure event
      if events.network[1][event[1]]
        # External exposure
        parentnode = 1
        branch_length = eventtimes[eventorder[i]]
      else
        # Internal exposure
        source = findfirst(events.network[2][:, event[1]])
        priorexposures = events.network[2][source, :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]])
        if isnan(eventtimes[source, 2]) || eventtimes[source, 2] > eventtimes[event[1], 1]
          # Undetected exposure source
          if any(priorexposures)
            # Prior exposures from this source
            parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
            branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
          else
            # No prior exposures from this source
            parentnode = get(eventnodes[source, 1])
            branch_length = eventtimes[eventorder[i]] - eventtimes[source, 1]
          end
        else
          # Detected exposure source
          if !any(priorexposures) || all(eventtimes[source, 2] .> eventtimes[priorexposures, 1])
            # Detection of exposure source is most recent, relevant event
            parentnode = get(eventnodes[source, 2])
            branch_length = eventtimes[eventorder[i]] - eventtimes[source, 2]
          else
            # Other, prior exposure is most recent, relevant event
            parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
            branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
          end
        end
      end
      branch!(tree, parentnode, branch_length)
      eventnodes[eventorder[i]] = length(tree.nodes)
    elseif event[2] == 2
      # Detection event
      priorexposures = events.network[2][event[1], :][:] & (eventtimes[:, 1] .< eventtimes[eventorder[i]])
      if any(priorexposures)
        # Individual has exposed others before detection
        parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
        branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
      else
        # Individual has not exposed others before detection
        parentnode = get(eventnodes[event[1], 1])
        branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
      end
      branch!(tree, parentnode, branch_length)
      eventnodes[eventorder[i]] = length(tree.nodes)
    elseif event[2] == 3
      # Removal event
      priorexposures = events.network[2][event[1], :][:]
      if any(priorexposures)
        # Individual has exposed others prior to removal
        if !isnan(eventtimes[event[1], 2]) && all(eventtimes[priorexposures, 1] .< eventtimes[event[1], 2])
          # Detection is most recent and relevant event
          parentnode = get(eventnodes[event[1], 2])
          branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 2]
        else
          # Exposure is the most recent and relevant event
          parentnode = get(eventnodes[priorexposures, 1][indmax(eventtimes[priorexposures, 1])])
          branch_length = eventtimes[eventorder[i]] - maximum(eventtimes[priorexposures, 1])
        end
      else
        # No prior exposures
        if !isnan(eventtimes[event[1], 2])
          # Has been previously detected
          parentnode = get(eventnodes[event[1], 2])
          branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 2]
        else
          # Has not been detected
          parentnode = get(eventnodes[event[1], 1])
          branch_length = eventtimes[eventorder[i]] - eventtimes[event[1], 1]
        end
      end
      branch!(tree, parentnode, branch_length)
      eventnodes[eventorder[i]] = length(tree.nodes)
    end
  end
  return tree, eventnodes[:,2]
end
