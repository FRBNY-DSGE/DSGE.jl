using Base.Test
using DSGE

include(joinpath(dirname(@__FILE__()),"util.jl"))

my_tests = [
            "core",
            "parameters",
            "solve/gensys",
            "solve/solve",
            "estimate/kalman",
            "estimate/posterior",
            "estimate/hessizero",
            "estimate/hessian",
            "estimate/csminwel",
            "estimate/optimize",
            "estimate/eig",
            "estimate/metropolis_hastings",
            "models/m990/m990",
            "models/smets_wouters/smets_wouters",
	    "data/data",
            # "forecast/smoothers",
            "forecast/filter",
            "forecast/forecast"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
