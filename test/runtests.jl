using Base.Test
using DSGE

include(joinpath(dirname(@__FILE__()),"util.jl"))

my_tests = [            
            ## "core",
            ## "models/m990/m990",
            ## "solve/gensys",
            ## "solve/solve",
            ## "solve/ordered_qz",
            ## "estimate/kalman",
            ## "estimate/posterior",
            ## #"estimate/hessian" (currently takes 6 hours to complete)
            ## "estimate/csminwel",
            ## "estimate/eig",
            "estimate/metropolis_hastings"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
