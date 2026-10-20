struct ColumnSusceptibility{V}
  values::V
end
(f::ColumnSusceptibility)(params, pop, i) = params[1] * f.values[i]

specialized_susceptibility(rf, rp, pop, i) = rf.susceptibility(rp.susceptibility, pop, i)
specialized_infectivity(rf, rp, pop, i, j) = rf.infectivity(rp.infectivity, pop, i, j)

@testset "Concrete risk storage and callable specialization" begin
  pop = Population(DataFrame(riskfactor1=[1.0, 2.0]), [0.0 3.0; 3.0 0.0])
  susceptibility = ColumnSusceptibility(pop.risks[!, :riskfactor1])
  for S in (SI, SIR, SEI, SEIR)
    funcs = Any[perf_constant, susceptibility, perf_distance_pair, perf_one]
    params = [[0.1], [2.0], Float64[], Float64[]]
    for state in (State_E, State_R)
      if state ∈ S
        push!(funcs, perf_constant)
        push!(params, [0.5])
      end
    end
    rf, rp = RiskFunctions{S}(funcs...), RiskParameters{S}(params...)
    @test all(isconcretetype, fieldtypes(typeof(rf)))
    @test all(isconcretetype, fieldtypes(typeof(rp)))
    @test (@inferred specialized_susceptibility(rf, rp, pop, 2)) == 4.0
    @test (@inferred specialized_infectivity(rf, rp, pop, 2, 1)) == 3.0
    sim = Simulation(pop, [State_I, State_S], rf, rp)
    @test sim isa Simulation
    @test sim.transmission_rates.internal[1, 2] == 12.0
  end
  @test all(isconcretetype, fieldtypes(typeof(pop)))
  @test Population(pop.risks).distances === nothing
  values = [0.2, 0.4]
  rp = RiskParameters{SI}(@view(values[1:1]), [1.0], Float64[], Float64[])
  @test rp.sparks isa SubArray
  @test rp.sparks[1] == 0.2
end
