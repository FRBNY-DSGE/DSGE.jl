using DSGE, JLD2, Distributions, PDMats, DataStructures, OrderedCollections, FileIO, Test, DataFrames, Dates, Nullables, Plots, Distributed, DelimitedFiles, Random

my_tests = [
#           "models/heterogeneous_agent/het_dsge_simple_taylor/het_dsge_simple_taylor",
           "models/heterogeneous_agent/het_dsge/het_dsge",
           #"models/heterogeneous_agent/het_dsge_lag/het_dsge_lag",
           "core",
           "parameters",
           #"models/heterogeneous_agent/het_dsge/het_dsge",
           #="models/representative_agent/an_schorfheide/an_schorfheide",
           "models/representative_agent/smets_wouters/smets_wouters",
           "models/representative_agent/m990/m990",
           "models/representative_agent/m1002/m1002",
           "models/representative_agent/m1010/m1010",
           "models/heterogeneous_agent/krusell_smith/krusell_smith",
           "models/heterogeneous_agent/bond_labor/bond_labor",=#
           "models/heterogeneous_agent/real_bond/real_bond",
           "models/heterogeneous_agent/real_bond_mkup/real_bond_mkup",
           # "models/heterogeneous_agent/krusell_smith_ct/krusell_smith_ct",
           # "models/heterogeneous_agent/one_asset_hank/one_asset_hank",
           # "models/heterogeneous_agent/one_asset_hank/interns",
           "data/misc",
           "data/load_data",
           "solve/gensys",
           "solve/solve",
           # "solve/solve_ct",
           # "solve/gensys_ct",
           # "solve/reduction",
           "estimate/filter",
           "estimate/cat",
           "estimate/posterior",
           "estimate/hessizero",
           "estimate/hessian",
           "estimate/csminwel",
           "estimate/optimize",
           "estimate/metropolis_hastings",
           # "estimate/smc/smc",
           "estimate/smc/helpers",
           "estimate/smc/initialization",
           "estimate/smc/resample",
           "estimate/smc/util",
           "estimate/smc/mutation",
           "forecast/smooth",
           "forecast/forecast",
           "forecast/shock_decompositions",
           "forecast/impulse_responses",
           "forecast/io",
           "forecast/forecast_one",
           "analysis/compute_meansbands",
           "altpolicy/altpolicy",
           "scenarios/scenario",
           "scenarios/forecast",
           "scenarios/switching",
           "scenarios/drivers",
           "decomp/decompose_forecast",
           "plot/plot"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
