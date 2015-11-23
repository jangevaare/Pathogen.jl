"""
Initial augmented data proposal
* For SEIR models, with infectivity rate, `ρ`
* With detection rate, `ν`
* Proposals made with exponential distributions
"""
function propose_augment(ρ::Float64,
                         ν::Float64,
                         obs::SEIR_observed,
                         debug=false::Bool)
  exposed_augmented = fill(NaN, length(obs.infectious))
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))
  for i = 1:length(obs.infectious)
    if !isnan(obs.infectious[i])
      infectious_augmented[i] = obs.infectious[i] - rand(Exponential(1/ν))
      exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
      if !isnan(obs.removed[i])
        removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0., obs.removed[i] - obs.infectious[i]))
      end
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Initial augmented data proposal
* For SEIR models, with infectivity rate, `ρ`
* With detection rate, `ν`
* Proposals made with uniform distributions
"""
function propose_augment(obs::SEIR_observed,
                         debug=false::Bool)
  exposed_augmented = fill(NaN, length(obs.infectious))
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))
  for i = 1:length(obs.infectious)
    if !isnan(obs.infectious[i])
      infectious_augmented[i] = rand(Uniform(0, obs.infectious[i]))
      exposed_augmented[i] = rand(Uniform(0, infectious_augmented[i]))
      if !isnan(obs.removed[i])
        removed_augmented[i] = rand(Uniform(obs.infectious[i], obs.removed[i]))
      end
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Initial augmented data proposal
* For SEIR models, with infectivity rate, `ρ`
* Proposals made with exponential distributions
"""
function propose_augment(ρ::Float64,
                         obs::SEIR_observed,
                         debug=false::Bool)
  infectious_augmented = obs.infectious
  removed_augmented = obs.removed
  exposed_augmented = fill(NaN, length(infectious_augmented))
  for i = 1:length(obs.infectious)
    if !isnan(infectious_augmented[i])
      exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Initial augmented data proposal
* For SEIR models, with infectivity rate, `ρ`
* Proposals made with uniform distributions
"""
function propose_augment(obs::SEIR_observed, debug=false::Bool)
  infectious_augmented = obs.infectious
  removed_augmented = obs.removed
  exposed_augmented = fill(NaN, length(infectious_augmented))
  for i = 1:length(obs.infectious)
    if !isnan(infectious_augmented[i])
      exposed_augmented[i] = rand(Uniform(0, infectious_augmented[i]))
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SEIR models with detection lag
* Requires `network` information
* Uses previous augmented data
* Proposals made with uniform distributions
* Generates new event times for a specified `changed_individuals`
"""
function propose_augment(changed_individuals::Vector{Int64},
                         network::Array{Bool, 2},
                         previous_aug::SEIR_augmented,
                         obs::SEIR_observed,
                         debug=false::Bool)
  exposed_augmented = previous_aug.exposed
  infectious_augmented = previous_aug.infectious
  removed_augmented = previous_aug.removed
  for i in changed_individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    pathway_in = pathwayto(i, network, debug)
    if debug
      if length(pathway_in) > 2
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection times (pathway from $i): $(infectious_augmented[pathway_out])")
        println("Augmented exposure times (pathway from $i): $(exposed_augmented[pathway_out])")
        println("Augmented infection time (exposer of $i): $(infectious_augmented[pathway_in[2]])")
        println("Augmented removal time (exposer of $i): $(removed_augmented[pathway_in[2]])")
      else
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection times (pathway from $i): $(infectious_augmented[pathway_out])")
        println("Augmented exposure times (pathway from $i): $(exposed_augmented[pathway_out])")
      end
    end
    # Randomize augmentation order
    aug_order = sample([1,2,3], 3, replace=false)
    for j in aug_order
      # Exposure time augmentation
      if j == 1
        if length(pathway_in) > 2
          if isnan(obs.removed[pathway_in[2]])
            exposed_augmented[i] = rand(Uniform(infectious_augmented[pathway_in[2]], infectious_augmented[i]))
          else
            exposed_augmented[i] = rand(Uniform(infectious_augmented[pathway_in[2]], minimum([infectious_augmented[i]; removed_augmented[pathway_in[2]]])))
          end
        else
          exposed_augmented[i] = rand(Uniform(0., infectious_augmented[i]))
        end
      # Infection time augmentation
      elseif j == 2
        infectious_augmented[i] = rand(Uniform(exposed_augmented[i], minimum([obs.infectious[i]; exposed_augmented[pathway_out[2:end]]])))
      # Removal time augmentation
      elseif !isnan(obs.removed[i])
        removed_augmented[i] = rand(Uniform(maximum([obs.infectious[i]; exposed_augmented[pathway_out[2:end]]]), obs.removed[i]))
      end
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SEIR models, with infectivity rate, `ρ`
* With detection rate, `ν`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified amount of `changes`
"""
function propose_augment(network::Array{Bool, 2},
                         previous_aug::SEIR_augmented,
                         obs::SEIR_observed,
                         changes=1::Int64,
                         debug=false::Bool)
  exposures = pathwayfrom(0, network, debug)[2:end]
  if changes == 0 || changes > length(exposures)
    changed_individuals = sample(exposures, length(exposures), replace=false)
  else
    changed_individuals = sample(exposures, changes, replace=false)
  end
  return propose_augment(changed_individuals, network, previous_aug, obs, debug)
end


"""
Proposes new augmented data
* For SEIR models, with infectivity rate, `ρ`
* With detection rate, `ν`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified `changed_individuals`
"""
function propose_augment(changed_individuals::Vector{Int64},
                         ρ::Float64,
                         ν::Float64,
                         network::Array{Bool, 2},
                         previous_aug::SEIR_augmented,
                         obs::SEIR_observed,
                         debug=false::Bool)
  exposed_augmented = previous_aug.exposed
  infectious_augmented = previous_aug.infectious
  removed_augmented = previous_aug.removed
  for i in changed_individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    pathway_in = pathwayto(i, network, debug)
    if debug
      if length(pathway_in) > 2
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection times (pathway from $i): $(infectious_augmented[pathway_out])")
        println("Augmented exposure times (pathway from $i): $(exposed_augmented[pathway_out])")
        println("Augmented infection time (exposer of $i): $(infectious_augmented[pathway_in[2]])")
        println("Augmented removal time (exposer of $i): $(removed_augmented[pathway_in[2]])")
      else
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection times (pathway from $i): $(infectious_augmented[pathway_out])")
        println("Augmented exposure times (pathway from $i): $(exposed_augmented[pathway_out])")
      end
    end
    # Randomize augmentation order
    aug_order = sample([1,2,3], 3, replace=false)
    for j in aug_order
      # Exposure time augmentation
      if j == 1
        if length(pathway_in) > 2
          if isnan(obs.removed[pathway_in[2]])
            exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), 0, infectious_augmented[i] - infectious_augmented[pathway_in[2]]))
          else
            exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), infectious_augmented[i] - minimum([infectious_augmented[i]; removed_augmented[pathway_in[2]]]), infectious_augmented[i] - infectious_augmented[pathway_in[2]]))
          end
        else
          exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
        end
      # Infection time augmentation
      elseif j == 2
        infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; exposed_augmented[pathway_out[2:end]]]), obs.infectious[i] - exposed_augmented[i]))
      # Removal time augmentation
      elseif !isnan(obs.removed[i])
        removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0, obs.removed[i] - maximum([obs.infectious[i]; exposed_augmented[pathway_out[2:end]]])))
      end
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SEIR models, with infectivity rate, `ρ`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified `changed_individuals`
"""
function propose_augment(changed_individuals::Vector{Int64},
                         ρ::Float64, network::Array{Bool, 2},
                         previous_aug::SEIR_augmented,
                         obs::SEIR_observed,
                         debug=false::Bool)
  exposed_augmented = previous_aug.exposed
  infectious_augmented = previous_aug.infectious
  removed_augmented = previous_aug.removed
  for i in changed_individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    pathway_in = pathwayto(i, network, debug)
    if debug
      if length(pathway_in) > 2
        println("Observed infection times (pathway from $i): $(infectious_augmented[pathway_out])")
        println("Augmented exposure times (pathway from $i): $(exposed_augmented[pathway_out])")
        println("Observed infection time (exposer of $i): $(infectious_augmented[pathway_in[2]])")
        println("Observed removal time (exposer of $i): $(removed_augmented[pathway_in[2]])")
      else
        println("Observed infection times (pathway from $i): $(infectious_augmented[pathway_out])")
        println("Augmented exposure times (pathway from $i): $(exposed_augmented[pathway_out])")
      end
    end
    if length(pathway_in) > 2
      if isnan(obs.removed[pathway_in[2]])
        exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), 0, infectious_augmented[i] - infectious_augmented[pathway_in[2]]))
      else
        exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), infectious_augmented[i] - minimum([infectious_augmented[i]; removed_augmented[pathway_in[2]]]), infectious_augmented[i] - infectious_augmented[pathway_in[2]]))
      end
    else
      exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Proposes augmented data
* For SEIR models, with infectivity rate, `ρ`
* With detection rate, `ν`
* Requires `network` information
* Proposals made with truncated exponential distributions
"""
function propose_augment(ρ::Float64,
                         ν::Float64,
                         network::Array{Bool, 2},
                         obs::SEIR_observed,
                         debug=false::Bool)
  individuals = pathwayfrom(0, network, debug)[2:end]
  exposed_augmented = fill(NaN, length(obs.infectious))
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))
  for i in individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    pathway_in = pathwayto(i, network, debug)
    if length(pathway_in) > 2
      if debug
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
        println("Augmented infection time (exposer of $i): $(infectious_augmented[pathway_in[2]])")
      end
      infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum(obs.infectious[pathway_out]), obs.infectious[i] - infectious_augmented[pathway_in[2]]))
      if isnan(obs.removed[pathway_in[2]])
        exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), 0., infectious_augmented[i] - infectious_augmented[pathway_in[2]]))
      else
        exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), infectious_augmented[i] - minimum([infectious_augmented[i]; removed_augmented[pathway_in[2]]]), infectious_augmented[i] - infectious_augmented[pathway_in[2]]))
      end
    else
      if debug
        println("Observed infection times (pathway from $i): $(obs.infectious[pathway_out])")
      end
      infectious_augmented[i] = obs.infectious[i] - rand(Truncated(Exponential(1/ν), obs.infectious[i] - minimum(obs.infectious[pathway_out]), Inf))
      exposed_augmented[i] = infectious_augmented[i] - rand(Exponential(1/ρ))
    end
    if !isnan(obs.removed[i])
      removed_augmented[i] = obs.removed[i] - rand(Truncated(Exponential(1/ν), 0., obs.removed[i] - obs.infectious[i]))
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SEIR models, with infectivity rate, `ρ`
* With detection rate, `ν`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified amount of `changes`
"""
function propose_augment(ρ::Float64,
                         ν::Float64,
                         network::Array{Bool, 2},
                         previous_aug::SEIR_augmented,
                         obs::SEIR_observed,
                         changes=1::Int64,
                         debug=false::Bool)
  exposures = pathwayfrom(0, network, debug)[2:end]
  if changes == 0 || changes > length(exposures)
    changed_individuals = sample(exposures, length(exposures), replace=false)
  else
    changed_individuals = sample(exposures, changes, replace=false)
  end
  return propose_augment(changed_individuals, ρ, ν, network, previous_aug, obs, debug)
end


"""
Proposes augmented data
* For SEIR models, with infectivity rate, `ρ`
* Requires `network` information
* Proposals made with truncated exponential distributions
"""
function propose_augment(ρ::Float64,
                         network::Array{Bool, 2},
                         obs::SEIR_observed,
                         debug=false::Bool)
  individuals = pathwayfrom(0, network, debug)[2:end]
  exposed_augmented = fill(NaN, length(obs.infectious))
  infectious_augmented = fill(NaN, length(obs.infectious))
  removed_augmented = fill(NaN, length(obs.removed))
  for i in individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    pathway_in = pathwayto(i, network, debug)
    infectious_augmented[i] = obs.infectious[i]
    if isnan(obs.removed[pathway_in[2]])
      exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), 0., infectious_augmented[i]-infectious_augmented[pathway_in[2]]))
    else
      exposed_augmented[i] = infectious_augmented[i] - rand(Truncated(Exponential(1/ρ), infectious_augmented[i] - minimum([infectious_augmented[i]; removed_augmented[pathway_in[2]]]), infectious_augmented[i] - infectious_augmented[pathway_in[2]]))
    end
    if !isnan(obs.removed[i])
      removed_augmented[i] = obs.removed[i]
    end
  end
  return SEIR_augmented(exposed_augmented, infectious_augmented, removed_augmented)
end


"""
Proposes new augmented data
* For SEIR models, with infectivity rate, `ρ`
* Requires `network` information
* Uses previous augmented data
* Proposals made with truncated exponential distributions
* Generates new event times for a specified amount of `changes`
"""
function propose_augment(ρ::Float64,
                         network::Array{Bool, 2},
                         previous_aug::SEIR_augmented,
                         obs::SEIR_observed,
                         changes=1::Int64,
                         debug=false::Bool)
  exposures = pathwayfrom(0, network, debug)[2:end]
  if changes == 0 || changes > length(exposures)
    changed_individuals = sample(exposures, length(exposures), replace=false)
  else
    changed_individuals = sample(exposures, changes, replace=false)
  end
  return propose_augment(changed_individuals, ρ, network, previous_aug, obs, debug)
end


"""
Log likelihood for detection rate, `ν`
* Requires `network` information
"""
function detection_loglikelihood(ν::Float64,
                                 aug::SEIR_augmented,
                                 network::Array{Bool, 2},
                                 obs::SEIR_observed,
                                 debug=false::Bool)
  ll = 0.
  individuals = find(!isnan(obs.infectious))
  exposed_augmented = aug.exposed
  infectious_augmented = aug.infectious
  removed_augmented = aug.removed
  for i in individuals
    pathway_out = pathwayfrom(i, network, 1, debug)
    ll += logpdf(Truncated(Exponential(1/ν), obs.infectious[i] - minimum([obs.infectious[i]; exposed_augmented[pathway_out[2:end]]]), obs.infectious[i] - exposed_augmented[i]), obs.infectious[i] - infectious_augmented[i])
    if !isnan(obs.removed[i])
      ll += logpdf(Truncated(Exponential(1/ν), 0, obs.removed[i] - maximum([obs.infectious[i]; exposed_augmented[pathway_out[2:end]]])), obs.removed[i] - removed_augmented[i])
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
