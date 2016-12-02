Busing Base.Test
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
            "estimate/smc",
            "estimate/mutation_RWMH",
            "models/m990/m990",
            "models/m1002/m1002",
            "models/AnSchorfheide/AnSchorfheide",
            "models/smets_wouters/smets_wouters",
            "data/misc",
            "data/load_data"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
