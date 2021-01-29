using ModelConstructors, Nullables, SMC, Test, Distributed, Distributions
using Dates, DataFrames, OrderedCollections, FileIO, DataStructures, LinearAlgebra, SparseArrays
using StatsBase, Random, CSV, StateSpaceRoutines, HDF5, JLD2, MAT, Plots
import ModelConstructors: @test_matrix_approx_eq, @test_matrix_approx_eq_eps
@everywhere using DSGE, JLD2, Printf, LinearAlgebra, ModelConstructors, SMC
HETDSGEGOVDEBT = "../src/models/heterogeneous/het_dsge_gov_debt/reference"

my_tests = [
            "abstractdsgemodel",
            "abstractvarmodel",
            "defaults",
            "parameters",
            "util",
            "statespace",

            "data/fred_data",
            "data/load_data",
            "data/load_data_poolmodel",
            "data/misc",
            "data/reverse_transform",
            "data/simulate_data",
            "data/transformations",
            "data/transform_data",
            "data/util",

            "solve/gensys",
            "solve/solve",
            "solve/gensys_uncertain_altpol",
            "solve/gensys2_uncertain_altpol_test1",
            "solve/gensys2_uncertain_altpol_test2",
            "solve/solve_poolmodel",
            #"solve/solve_ct",
            #"solve/gensys_ct",
            #"solve/reduction",

            "estimate/filter",
            "estimate/cat",
            "estimate/posterior",
            "estimate/poolmodel_tpf",
            "estimate/filter_poolmodel",
            "estimate/posterior_poolmodel",
            "estimate/estimate_bma",
            "estimate/hessian",
            "estimate/util",
            "estimate/csminwel",
            "estimate/optimize",
            "estimate/var/dsgevar_likelihood",
            "estimate/var/dsgevecm_likelihood",

            "estimate/metropolis_hastings",
            "estimate/regime_switching_mh",
            "estimate/smc/helpers",
            "estimate/smc/initialization",
            "estimate/smc/mutation",
            "estimate/smc/resample",
            "estimate/smc/particle",
            "estimate/smc/smc",
            "estimate/smc/regime_switching_smc",

            "forecast/drivers",
            "forecast/smooth",
            "forecast/forecast",
            "forecast/shock_decompositions",
            "forecast/impulse_responses",
            "forecast/io",
            "forecast/forecast_one",
            "forecast/automatic_tempalt_zlb",
            "forecast/forecast_regime_switching",
            "forecast/time_varying_credibility",
            "forecast/util",
            "forecast/var/impulse_responses",
            "forecast/var/dsgevar/impulse_responses",
            "forecast/var/dsgevecm/impulse_responses",

            "analysis/compute_meansbands",
            "analysis/df_to_table",
            "analysis/io",
            "analysis/meansbands",
            "analysis/moments",
            "analysis/util",

            "altpolicy/altpolicy",

            "scenarios/scenario",
            "scenarios/forecast",
            "scenarios/switching",
            "scenarios/drivers",

            "decomp/decompose_forecast",

            "plot/util",

            "models/representative/smets_wouters/smets_wouters",
            "models/representative/smets_wouters_orig/smets_wouters_orig",
            "models/representative/m990/m990",
            "models/representative/m1002/m1002",
            "models/representative/m1010/m1010",
            "models/representative/m904/m904",
            "models/representative/m805/m805",
            "models/poolmodel/poolmodel",
            "models/var/dsgevar/dsgevar",
            "models/var/dsgevecm/dsgevecm",
            "models/var/util",
            "models/heterogeneous/het_dsge_gov_debt/het_dsge_gov_debt_reduce_ell"
            # "models/heterogeneous/het_dsge_gov_debt/het_dsge_gov_debt"
            # "models/representative/rep_dsge_gov_debt/rep_dsge_gov_debt",
            # "models/heterogeneous/het_dsge_simple_taylor/het_dsge_simple_taylor",
            # "models/heterogeneous/het_dsge/het_dsge",
            # "models/heterogeneous/het_dsge_lag/het_dsge_lag",
            # "models/heterogeneous/het_dsge/het_dsge",

            # "models/heterogeneous/krusell_smith/krusell_smith",
            # "models/heterogeneous/bond_labor/bond_labor",
            # "models/heterogeneous/real_bond/real_bond",
            # "models/heterogeneous/real_bond_mkup/real_bond_mkup",
            # "models/heterogeneous/krusell_smith_ct/krusell_smith_ct",
            # "models/heterogeneous/one_asset_hank/one_asset_hank",
            # "models/heterogeneous/one_asset_hank/interns",
            ]

if VERSION >= v"1.3"
    my_tests = vcat([
                     "packet/packet", # These two tests generate segmentation fault errors
                     "plot/plot",     # in lower versions of Julia (at least w/1.0 and 1.1)
                    ],
                    my_tests)
end

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
