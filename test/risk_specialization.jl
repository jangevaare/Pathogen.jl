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
    samples = [rp, rp]
    expected = repeat(permutedims(convert(Vector{Float64}, rp)), 2, 1)
    @test convert(Matrix{Float64}, samples) == expected
    recipe = Pathogen.RecipesBase.apply_recipe(Dict{Symbol, Any}(), 1:2, samples)
    @test recipe[1].args[2] == expected
  end
  @test all(isconcretetype, fieldtypes(typeof(pop)))
  @test Population(pop.risks).distances === nothing
  values = [0.2, 0.4]
  rp = RiskParameters{SI}(@view(values[1:1]), [1.0], Float64[], Float64[])
  @test rp.sparks isa SubArray
  @test rp.sparks[1] == 0.2
  parameter_free = RiskParameters{SI}(nothing, nothing, nothing, nothing)
  rf = RiskFunctions{SI}(perf_one, perf_one, (p, pop, i, j) -> 0.0, perf_one)
  sim = Simulation(pop, [State_I, State_S], rf, parameter_free)
  @test sim.event_rates.infection == [0.0, 1.0]
end
