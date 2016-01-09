"""
Initial augmented data proposal
* For SIR models
* With detection rate, `ν`
* Proposals made with exponential distributions
"""
function propose_augment(ν::Float64, obs::SIR_observed, debug=false::Bool)
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))
  for i = 1:length(obs.infectious)
    if !isnan(obs.infectious[i])
      infectious_augmented[i] = obs.infectious[i] - rand(Exponential(1/ν))
      if !isnan(obs.removed[i])
        removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0., obs.removed[i] - obs.infectious[i]))
      end
    end
  end
  return SIR_augmented(infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SIR models
* With detection rate, `ν`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified `changed_individuals`
"""
function propose_augment(changed_individuals::Vector{Int64},
                         ν::Float64,
                         network::Array{Bool, 2},
                         previous_aug::SIR_augmented,
                         obs::SIR_observed,
                         debug=false::Bool)
  infectious_augmented = previous_aug.infectious
  removed_augmented = previous_aug.removed
  for i in changed_individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    pathway_in = pathwayto(i, network, debug)
    if debug
      if length(pathway_in) > 2
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection times (pathway from $i): $(infectious_augmented[pathway_out])")
        println("Augmented infection time (exposer of $i): $(infectious_augmented[pathway_in[2]])")
        println("Augmented removal time (exposer of $i): $(removed_augmented[pathway_in[2]])")
      else
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection times (pathway from $i): $(infectious_augmented[pathway_out])")
      end
    end
    # Randomize augmentation order
    aug_order = sample([1,2], 2, replace=false)
    for j in aug_order
      # Exposure time augmentation
      if j == 1
        if length(pathway_in) > 2
          if isnan(obs.removed[pathway_in[2]])
            infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]]), obs.infectious[i] - infectious_augmented[pathway_in[2]]))
          else
            infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]; removed_augmented[pathway_in[2]]]), obs.infectious[i] - infectious_augmented[pathway_in[2]]))
          end
        else
          infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]]), Inf))
        end
      # Removal time augmentation
      elseif j==2 && !isnan(obs.removed[i])
        removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0., obs.removed[i] - maximum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]])))
      end
    end
  end
  return SIR_augmented(infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SIR models
* With detection rate, `ν`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified `changed_individual`
"""
function propose_augment(i::Int64,
                         ν::Float64,
                         network::Array{Bool, 2},
                         previous_aug::SIR_augmented,
                         obs::SIR_observed,
                         debug=false::Bool)
  infectious_augmented = previous_aug.infectious
  removed_augmented = previous_aug.removed
  pathway_out = pathwayfrom(i, network, 1, debug)
  pathway_in = pathwayto(i, network, debug)
  if debug
    println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
    println("Augmented infection times (pathway from $i): $(infectious_augmented[pathway_out])")
  end
  # Exposure time augmentation
  infectious_augmented[i] = obs.infectious[i] - rand(Exponential(1/ν))
  difference = infectious_augmented[i] - previous_aug.infectious[i]
  for j in pathway_out[2:end]
    infectious_augmented[j] += difference
    removed_augmented[j] += difference
  end

  # Removal time augmentation
  if !isnan(obs.removed[i])
    removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0., obs.removed[i] - maximum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]])))
  end
  return SIR_augmented(infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SIR models
* With detection rate, `ν`
* Requires `network` information
* Proposals made with truncated exponential distributions
"""
function propose_augment(ν::Float64,
                         network::Array{Bool, 2},
                         obs::SIR_observed,
                         debug=false::Bool)

  infectious_augmented = obs.infectious
  removed_augmented = obs.removed
  exposures = pathwayfrom(0, network, debug)[2:end]
  for i in exposures
    pathway_out = pathwayfrom(i, network, 1, debug)
    pathway_in = pathwayto(i, network, debug)

    if debug
      if length(pathway_in) > 2
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection time (exposer of $i): $(infectious_augmented[pathway_in[2]])")
        println("Augmented removal time (exposer of $i): $(removed_augmented[pathway_in[2]])")
      else
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
      end
    end

    # Infection time augmentation
    if length(pathway_in) > 2
      if isnan(obs.removed[pathway_in[2]])
        infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]]), obs.infectious[i] - infectious_augmented[pathway_in[2]]))
      else
        infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]; removed_augmented[pathway_in[2]]]), obs.infectious[i] - infectious_augmented[pathway_in[2]]))
      end
    else
      infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]]), Inf))
    end

    # Removal time augmentation
    if !isnan(obs.removed[i])
      removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0., obs.removed[i] - maximum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]])))
    end

  end
  return SIR_augmented(infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SIR models
* With detection rate, `ν`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified amount of `changes`
"""
function propose_augment(ν::Float64,
                         network::Array{Bool, 2},
                         previous_aug::SIR_augmented,
                         obs::SIR_observed,
                         changes=rand(Poisson(1.))::Int64,
                         debug=false::Bool)
  exposures = pathwayfrom(0, network, debug)[2:end]
  if changes == 0 || changes > length(exposures)
    return propose_augment(ν, obs, debug)
  else
    changed_individuals = sample(exposures, changes, replace=false)
    return propose_augment(changed_individuals, ν, network, previous_aug, obs, debug)
  end
end


"""
Log likelihood for detection rate, `ν` for SIR models
* Requires `network` information
"""
function detection_loglikelihood(ν::Float64,
                                 aug::SIR_augmented,
                                 network::Array{Bool, 2},
                                 obs::SIR_observed,
                                 debug=false::Bool)
  ll = 0.
  individuals = find(!isnan(obs.infectious))
  infectious_augmented = aug.infectious
  removed_augmented = aug.removed
  for i in individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    if length(pathway_in) > 2
      if isnan(obs.removed[pathway_in[2]])
        ll += logpdf(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]]), obs.infectious[i] - infectious_augmented[pathway_in[2]]), obs.infectious[i] - infectious_augmented[i])
      else
        ll += logpdf(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]; removed_augmented[pathway_in[2]]]), obs.infectious[i] - infectious_augmented[pathway_in[2]]), obs.infectious[i] - infectious_augmented[i])
      end
    else
      ll += logpdf(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]]), Inf), obs.infectious[i] - infectious_augmented[i])
    end
    if !isnan(obs.removed[i])
      ll += logpdf(Truncated(Exponential(1/ν), 0., obs.removed[i] - maximum([obs.infectious[i]; infectious_augmented[pathway_out[2:end]]])), obs.removed[i] - removed_augmented[i])
    end
    if isnan(ll)
      ll = -Inf
    end
    ll == -Inf && break
  end
  if debug
    println("Detection log likelihood: $(round(ll, 3))")
  end
  return ll
end


"""
Initial augmented data proposal
* For SIR models
"""
function propose_augment(obs::SIR_observed, debug=false::Bool)
  infectious_augmented = obs.infectious
  removed_augmented = obs.removed
  return SIR_augmented(infectious_augmented, removed_augmented)
end
