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

@time hessian, stoph = hessian!(model, mode, YY; verbose=true)
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
