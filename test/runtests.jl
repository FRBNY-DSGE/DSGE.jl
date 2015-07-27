using Base.Test
using DSGE

include("util.jl")

my_tests = [
            "abstractmodel",
            "models/m990/m990",
            "solve/solve",
            "solve/ordered_qz",
            "estimate/kalman",
            "estimate/posterior"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
