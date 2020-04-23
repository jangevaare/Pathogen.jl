module Pathogen

  # Dependencies
  using Distributed,
        DataFrames,
        Distributions,
        RecipesBase,
        Logging,
        StatsBase,
        Statistics,
        ProgressMeter,
        LinearAlgebra,
        OnlineStats

  import StatsBase.mean,
         StatsBase.mode,
         Base.summary,
         ProgressMeter.next!

  # Core
  include("core/EpidemicModel.jl")
  include("core/DiseaseState.jl")
  include("core/Population.jl")
  include("core/Events/Event.jl")
  include("core/Events/Events.jl")
  include("core/Events/EventRates.jl")
  include("core/Events/EventObservations.jl")
  include("core/Risks/AbstractRisk.jl")
  include("core/Risks/RiskFunctions.jl")
  include("core/Risks/RiskParameters.jl")
  include("core/Transmissions/Transmission.jl")
  include("core/Transmissions/AbstractTransmissionNetwork.jl")
  include("core/Transmissions/TransmissionNetwork.jl")
  include("core/Transmissions/TransmissionRates.jl")
  include("core/initialize.jl")
  include("core/update!.jl")

  # Simulation
  include("simulation/Simulation.jl")
  include("simulation/generate.jl")
  include("simulation/simulate!.jl")
  include("simulation/update!.jl")
  include("simulation/observe.jl")

  # Inference
  include("inference/TransmissionNetworkDistribution.jl")
  include("inference/RiskPriors.jl")
  include("inference/EventExtents.jl")
  include("inference/MarkovChain.jl")
  include("inference/MCMC.jl")
  include("inference/_pathway_from.jl")
  include("inference/_pathway_to.jl")
  include("inference/_accept.jl")
  include("inference/_bounds.jl")
  include("inference/generate.jl")
  include("inference/initialize.jl")
  include("inference/update!.jl")
  include("inference/loglikelihood.jl")
  include("inference/logpriors.jl")
  include("inference/start!.jl")
  include("inference/iterate!.jl")
  include("inference/summary.jl")

  # Visualization
  include("visualization/epidemic_curve.jl")
  include("visualization/observations.jl")
  include("visualization/population.jl")
  include("visualization/transmission_network.jl")
  include("visualization/transmission_network_distribution.jl")
  include("visualization/degree_distribution.jl")
  include("visualization/trace.jl")

  export
    #re-exports from DataFrames.jl
    DataFrame,
    #re-exports from Distributions.jl
    UnivariateDistribution,
    Arcsine,
    Bernoulli,
    Beta,
    BetaBinomial,
    BetaPrime,
    Binomial,
    Biweight,
    Categorical,
    Cauchy,
    Chernoff,
    Chi,
    Chisq,
    Cosine,
    DiagNormal,
    DiagNormalCanon,
    Dirichlet,
    DirichletMultinomial,
    DiscreteUniform,
    DoubleExponential,
    EdgeworthMean,
    EdgeworthSum,
    EdgeworthZ,
    Erlang,
    Epanechnikov,
    Exponential,
    FDist,
    FisherNoncentralHypergeometric,
    Frechet,
    FullNormal,
    FullNormalCanon,
    Gamma,
    DiscreteNonParametric,
    GeneralizedPareto,
    GeneralizedExtremeValue,
    Geometric,
    Gumbel,
    Hypergeometric,
    InverseWishart,
    InverseGamma,
    InverseGaussian,
    IsoNormal,
    IsoNormalCanon,
    Kolmogorov,
    KSDist,
    KSOneSided,
    Laplace,
    Levy,
    LKJ,
    LocationScale,
    Logistic,
    LogNormal,
    LogitNormal,
    MatrixBeta,
    MatrixFDist,
    MatrixNormal,
    MatrixTDist,
    MixtureModel,
    Multinomial,
    MultivariateNormal,
    MvLogNormal,
    MvNormal,
    MvNormalCanon,
    MvNormalKnownCov,
    MvTDist,
    NegativeBinomial,
    NoncentralBeta,
    NoncentralChisq,
    NoncentralF,
    NoncentralHypergeometric,
    NoncentralT,
    Normal,
    NormalCanon,
    NormalInverseGaussian,
    Pareto,
    PGeneralizedGaussian,
    Product,
    Poisson,
    PoissonBinomial,
    QQPair,
    Rayleigh,
    Semicircle,
    Skellam,
    StudentizedRange,
    SymTriangularDist,
    TDist,
    TriangularDist,
    Triweight,
    Truncated,
    Uniform,
    UnivariateGMM,
    VonMises,
    VonMisesFisher,
    WalleniusNoncentralHypergeometric,
    Weibull,
    Wishart,
    ZeroMeanIsoNormal,
    ZeroMeanIsoNormalCanon,
    ZeroMeanDiagNormal,
    ZeroMeanDiagNormalCanon,
    ZeroMeanFullNormal,
    ZeroMeanFullNormalCanon,
    # Pathogen.jl types/functions
    EpidemicModel,
    SEIR, SEI, SIR, SI,
    DiseaseState,
    State_S, State_E, State_I, State_R,
    RiskFunctions, RiskParameters, RiskPriors,
    Population,
    AbstractTransmissionNetwork,
    TransmissionNetwork,
    TransmissionNetworkDistribution, TNDistribution,
    TransmissionNetworkPrior, TNPrior,
    TransmissionNetworkPosterior, TNPosterior,
    EndogenousTransmission, ExogenousTransmission, NoTransmission,
    Simulation,
    next!, simulate!,
    Events, EventObservations, EventExtents,
    observe,
    MCMC, start!, iterate!,
    summary, mean, mode
end
