using Base.Test
using DSGE

include(joinpath(dirname(@__FILE__()),"util.jl"))

my_tests = [            
            ## "core",
            ## "parameters",
            ## "solve/gensys",
            ## "solve/solve",
            ## "solve/ordered_qz",
            ## "estimate/kalman",
            ## "estimate/posterior",
            ## "estimate/hessizero",
            ## "estimate/hessian",
            ## "estimate/csminwel",
            ## "estimate/optimize",
            ## "estimate/eig",
            ## "models/m990/m990",
            "estimate/metropolis_hastings"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
