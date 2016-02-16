using Base.Test
using DSGE
using Distributions

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
            "forecast/smoothers",
            "models/m990/m990",
            "models/smets_wouters/smets_wouters"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
