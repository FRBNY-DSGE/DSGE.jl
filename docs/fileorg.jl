# This file has extensoion .jl only for the purposes of syntax highlighting; it should not
#   be run

dsge.jl/
  src/
    DSGE.jl
      module DSGE
        include("init/DistributionsExt.jl")
        include("init/FinancialFrictionsFunctions.jl")
        include("solve/Gensys.jl")
        using Distributions, Roots, MATLAB
        using .DistributionsExt, .Gensys
        import Base: convert, promote_rule, log, exp, start, next, done
        export AbstractModel, Param, toreal, tomodel, Parameters, prior, ModelInds, makedict, solve, dlyap, kalcvf2NaN, kalsmth_k93, likelihood, Parameters990, Model990, model_specifications, eqcond, measurement
        include("core.jl")
        include("solve/solve.jl")
        include("estimate/kalman.jl")
        include("estimate/likelihood.jl")
        include("models/m990/m990.jl")
        include("models/m990/parameters.jl")
        include("models/m990/modelinds.jl")
      end
    core.jl
      abstract AbstractModel
      type Param
      function toreal
      function tomodel
      abstract Parameters
      function prior
      type ModelInds
      function makedict
    models/
      m990/
        m990.jl
          type Model990 <: AbstractModel
          function model_specifications
          function eqcond
          function measurement
        parameters.jl
          type Parameters990 <: Parameters
          function steadystate!
        modelinds.jl
          function ModelInds
    init/
      DistributionsExt.jl
        module DistributionsExt
          using Distributions
          import Distributions: pdf, mean, std
          type PointMass <: Distribution{Univariate, Continuous}
          function Beta
          function Gamma
          function InverseGamma
        end
      FinancialFrictionsFunctions.jl
        module FinancialFrictionsFunctions
          using Distributions
          export ζ_spb_fn, ...
          function ζ_spb_fn
          ...
        end
    solve/
      Gensys.jl
        module Gensys
          include("ordered_qz.jl")
          export gensys, ordschur
          function new_div
          function gensys
          function gensys!
          function qzdiv
          function qzdiv!
          function qzswitch!
        end
      Gensys_versions.jl
        function gensys_qzdiv
        function gensys_ordschur
      ordered_qz.jl
      solve.jl
        using ..AbstractModel, ..Gensys
        function solve
        function augment_states
    estimate/
      kalman.jl
        function kalcvf2NaN
        function kalsmth_k93
        function distsmth_k93
      likelihood.jl
        function likelihood
        function dlyap
    forecast/
  docs/
    fileorg.jl # You are here
    ProposedOrganization.md
    DSGE_Model_Documentation.pdf
    supersticky1214.pdf
  test/
    test_readcsv_complex.csv
    util.jl
      using Base.Test
      function minusnan
      function test_matrix_eq
      function readcsv_complex
      function test_util
      function test_test_matrix_eq
      function test_readcsv_complex
    test_AbstractModel.jl
      using Base.Test, Distributions
      using DSGE
      using DSGE: DistributionsExt
      function test_all
      function test_param
      function test_parameters
      function test_modelinds
    models/
      test_M990.jl
        using Base.Test, Distributions
        using DSGE
        using DSGE: DistributionsExt
        include("../util.jl")
        function test_all
        function test_parameters
        function test_modelinds
        function test_eqcond
        function test_measurement
      m990/
        parameters/
          para.csv
          ...
        eqcond/
          Γ0.csv
          ...
    solve/
      test_Gensys.jl
        using Base: Test, LinAlg
        using MATLAB
        using DSGE: Gensys
        include("../util.jl")
        include("../../src/solve/Gensys_versions.jl")
      test_ordered_qz.jl
        using Base.Test, Compat
        import Base.Linalg: BlasComplex, BlasFloat, BlasReal, QRPivoted
        include("../../src/solve/ordered_qz.jl")
      gensys/
        gensys_args.mat
        gensys_variables.mat
        make_gensys_args.m
        test_Gensys.m
      test_solve.mat
      test_solve.jl
        using Base: Test
        using MATLAB
        using DSGE: M990
        include("../util.jl")
    estimate/
      test_Kalman.jl
        using Base.Test
        using MATLAB
        using DSGE.Kalman
        include("../util.jl")
      kalcvf2NaN/
        kalcvf2NaN_args.mat
        kalcvf2NaN_out7.mat
        kalcvf2NaN_out9.mat
        test_kalcvf2NaN.m
      test_likelihood.mat
      test_likelihood.jl
        using Base: Test
        using MATLAB
        using DSGE
        using DSGE: M990
        include("../util.jl")
  README.md
