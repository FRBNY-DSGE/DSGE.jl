# This file has extension .jl only for the purposes of syntax highlighting; it should not be run

dsge.jl/
  src/
    DSGE.jl
      module DSGE
        export AbstractModel, DistributionsExt, FinancialFrictionsFunctions, M990
        include("init/DistributionsExt.jl")
        include("init/FinancialFrictionsFunctions.jl")
        include("AbstractModel.jl")
        include("models/M990.jl")
      end
    AbstractModel.jl
      module AbstractModel
        export Param, toreal, tomodel, Parameters, logprior, ModelInds, makedict, Model
        include("abstractmodel/param.jl")
        include("abstractmodel/parameters.jl")
        include("abstractmodel/modelinds.jl")
        include("abstractmodel/model.jl")
      end
    abstractmodel/
      param.jl
        using Distributions: Distribution
        import Base: convert, promote_rule, log, exp
        using DSGE.DistributionsExt: PointMass
        type Param
        function toreal
        function tomodel
      parameters.jl
        abstract Parameters
        function logprior
      modelinds.jl
        type ModelInds
        function makedict
      model.jl
        type Model
    models/
      M990.jl
        module M990
          using DSGE.AbstractModel
          export Parameters990, ModelInds990, Model990, eqcond
          include("spec.jl")
          include("parameters.jl")
          include("modelinds.jl")
          include("eqcond.jl")
          function Model990
        end
      m990/
        spec.jl
        parameters.jl
          using Distributions: Normal, quantile
          using Roots: fzero
          using DSGE.DistributionsExt: Beta, Gamma, InverseGamma
          using DSGE.FinancialFrictionsFunctions
          type Parameters990 <: Parameters
          function steadystate!
        modelinds.jl
          using DSGE.AbstractModel
          function ModelInds990
        eqcond.jl
          function eqcond
    init/
      DistributionsExt.jl
        module DistributionsExt
          using Distributions
          import Distributions: pdf
          type PointMass <: Distribution{Univariate, Continuous}
          function Distributions.pdf
          function Beta
          function Gamma
          function InverseGamma
        end
      FinancialFrictionsFunctions.jl
        module FinancialFrictionsFunctions
          using Distributions
          export zetaspbfcn, ...
          function zetaspbfcn
          ...
        end
    solve/
    estimate/
    forecast/
  docs/
    fileorg.jl # You are here
    ProposedOrganization.md
    DSGE_Model_Documentation.pdf
    supersticky1214.pdf
  test/
    util.jl
      function test_matrix_eq
      function test_util
    test_AbstractModel.jl
      using Base.Test
      using Distributions
      using DSGE: DistributionsExt, AbstractModel
      function test_all
      function test_param
      function test_parameters
      function test_modelinds
    models/
      test_M990.jl
        using Base.Test
        using DSGE: DistributionsExt, AbstractModel, M990
        include("../util.jl")
        function test_all
        function test_parameters
        function test_modelinds
        function test_eqcond
      m990/
        G0.csv
        G1.csv
        C.csv
        PSI.csv
        PIE.csv
  README.md
