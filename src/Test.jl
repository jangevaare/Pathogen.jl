sim_gene=rand(1:4,700,1)
time = [0.0,]
time_max = 1
i = 1
while time[end] < time_max
  mutation_rates=sum(JC69((0.5,)),1)[sim_gene[:,i]]
  mutation_time = rand(Exponential(1/sum(mutatation_rates)))
  if mutation_time + time[i] < time_max
    push!(time, time[i]+mutation_time)
    mutation = find(rand(Multinomial(1, mutatation_rates/sum(mutatation_rates))))[1]
    sim_gene = hcat(sim_gene, sim_gene[:,i])
    sim_gene[mutation, i+1] = find(rand(Multinomial(1, JC69((0.5,))[sim_gene[mutation],:][:]/mutation_rates[mutation])))[1]
    i += 1
  else
    push!(time, time_max)
    sim_gene = hcat(sim_gene, sim_gene[:,i])
    i += 1
  end
end

time

sim_gene
