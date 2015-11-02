using Base: Test
using DSGE, HDF5

include("../util.jl")
path = dirname(@__FILE__)

mode = h5open("$path/../reference/mode_in.h5") do file
    read(file, "params")
end
YY = h5open("$path/../reference/YY.h5") do file
    read(file, "YY")
end
hessian_expected = h5open("$path/../reference/hessian.h5") do file
    read(file, "hessian")
end

model = Model990()

@time hessian, hessian_errors = hessian!(model, mode, YY; verbose=true)
@test test_matrix_eq(hessian_expected, hessian; ϵ=1.0, verbose=true)

# Output 2015-10-30 - normal code
# 3893.685200 seconds (1.64 G allocations: 4.428 TB, 5.56% gc time)
# 734 of 7569 entries have opposite signs
# 300 of 7569 entries have different signs, but one of the entries is 0
# 4084 of 7569 entries with abs diff > 0
# 0 of 7569 entries with abs diff > 1.0
# Max abs diff of 0.5182971716154725 at entry (11,9)
# The entries at (11,9) are -360.19713132228446 + 0.0im (expected) and -359.678834150669 + 0.0im (actual)
# Max percent error of 56.59815003314423 at entry (68,3)
# The entries at (68,3) are 2.7578448419364482e-5 + 0.0im (expected) and -0.0015333107129011558 + 0.0im (actual)

# Output - hessiandiagelement function only used
# 3925.981645 seconds (1.64 G allocations: 4.428 TB, 5.59% gc time)
# 734 of 7569 entries have opposite signs
# 300 of 7569 entries have different signs, but one of the entries is 0
# 4084 of 7569 entries with abs diff > 0
# 0 of 7569 entries with abs diff > 1.0
# Max abs diff of 0.5182971716154725 at entry (11,9)
# The entries at (11,9) are -360.19713132228446 + 0.0im (expected) and -359.678834150669 + 0.0im (actual)
# Max percent error of 56.59815003314423 at entry (68,3)
# The entries at (68,3) are 2.7578448419364482e-5 + 0.0im (expected) and -0.0015333107129011558 + 0.0im (actual)

# Output - commit c723*, a lot of refactoring.
# 3913.768002 seconds (1.64 G allocations: 4.428 TB, 5.78% gc time)
# 734 of 7569 entries have opposite signs
# 300 of 7569 entries have different signs, but one of the entries is 0
# 4084 of 7569 entries with abs diff > 0
# 0 of 7569 entries with abs diff > 1.0
# Max abs diff of 0.5182971716154725 at entry (11,9)
# The entries at (11,9) are -360.19713132228446 + 0.0im (expected) and -359.678834150669 + 0.0im (actual)
# Max percent error of 56.59815003314423 at entry (68,3)
# The entries at (68,3) are 2.7578448419364482e-5 + 0.0im (expected) and -0.0015333107129011558 + 0.0im (actual)

# Output (15 workers) - commit f086, everything parallelized
# 337.695092 seconds (3.93 M allocations: 197.360 MB, 0.01% gc time)
# 46 of 7569 entries have opposite signs
# 36 of 7569 entries have different signs, but one of the entries is 0
# 3484 of 7569 entries with abs diff > 0
# 0 of 7569 entries with abs diff > 1.0
# Max abs diff of 0.6220439783200788 at entry (17,10)
# The entries at (17,10) are 37.70649227325377 + 0.0im (expected) and 38.32853625157385 + 0.0im (actual)
# Max percent error of 54.59815003314423 at entry (68,32)
# The entries at (68,32) are -2.7578448419364482e-5 + 0.0im (expected) and -0.0015333107129011558 + 0.0im (actual)
