using Base.Test
using DSGE

include("util.jl")

DSGE_TESTPATH = dirname(@__FILE__())
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
            "estimate/metropolis_hastings"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)

    ## test_path = joinpath(DSGE_TESTPATH, test_file)
    ## color = Base.have_color? "--color=yes" : "--color=no"
    ## #codecov = coverage? ["--code-coverage=user", "--inline=no"] : ["--code-coverage=none"]
    ## codecov =  ["--code-coverage=none"]
    ## julia_exe = joinpath(JULIA_HOME, Base.julia_exename())
    ## run(`$julia_exe --check-bounds=yes $codecov $color $test_path`)
end
