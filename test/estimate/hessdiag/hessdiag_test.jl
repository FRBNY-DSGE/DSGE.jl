using DSGE, MATLAB

ndx = 6
dx = exp(-(6:2:(6+(ndx-1)*2))')

# For each i, j = 1:82, k = 1:6, `hessizero` calculates
#     hessdiag[k] = -(fx - fdx - fdy + fdxdy)/(dx[k])^2
#     hessian[i, j] = -(hessdiag[3] + hessdiag[4])/2
# hessdiag_mat[i, j, k] and hessdiag_jl[i, j, k] correspond to hessdiag[k] for entry (i, j)

mf = MatFile("hessdiag.mat")
hessdiag_mat = get_variable(mf, "hessdiag_mat")
hessdiag_jl = get_variable(mf, "hessdiag_jl")
close(mf)

# Hypothesis: As dx[k] gets bigger, hessdiag_mat[i, j, k] and hessdiag_jl[i, j, k] get closer
# Hence hessian_12_mat and hessian_12_jl should be closer to each other than hessian_56_mat and
#   hessian_56_jl
# This is because the posterior distribution is relatively flat and small discrepancies between the
#   Matlab and Julia f-values get magnified when divided by (dx[k])^2

hessian_12_mat = [-(hessdiag_mat[i, j, 1] + hessdiag_mat[i, j, 2])/2 for i = 1:82, j = 1:82]
hessian_34_mat = [-(hessdiag_mat[i, j, 3] + hessdiag_mat[i, j, 4])/2 for i = 1:82, j = 1:82]
hessian_56_mat = [-(hessdiag_mat[i, j, 5] + hessdiag_mat[i, j, 6])/2 for i = 1:82, j = 1:82]
hessian_12_jl = [-(hessdiag_jl[i, j, 1] + hessdiag_jl[i, j, 2])/2 for i = 1:82, j = 1:82]
hessian_34_jl = [-(hessdiag_jl[i, j, 3] + hessdiag_jl[i, j, 4])/2 for i = 1:82, j = 1:82]
hessian_56_jl = [-(hessdiag_jl[i, j, 5] + hessdiag_jl[i, j, 6])/2 for i = 1:82, j = 1:82]

# Absolute differences
abs_diff(x, y) = abs(abs(x) - abs(y))

abs_diffs_12 = map(abs_diff, hessian_12_mat, hessian_12_jl)
abs_diffs_34 = map(abs_diff, hessian_34_mat, hessian_34_jl)
abs_diffs_56 = map(abs_diff, hessian_56_mat, hessian_56_jl)

mean(abs_diffs_12)
mean(abs_diffs_34)
mean(abs_diffs_56)

median(abs_diffs_12)
median(abs_diffs_34)
median(abs_diffs_56)

maximum(abs_diffs_12)
maximum(abs_diffs_34)
maximum(abs_diffs_56)

abs_dict_bool = [(i, j) => abs_diffs_12[i, j] <= abs_diffs_34[i, j] <= abs_diffs_56[i, j] for i = 1:82, j = 1:82]
count(x -> x, values(abs_dict_bool))

# Percent differences
function pct_diff(x, y)
    diff = abs_diff(x, y)
    if x != 0
        return diff/abs(x)
    elseif y != 0
        return diff/abs(y)
    else
        return 0
    end
end

pct_diffs_12 = map(pct_diff, hessian_12_mat, hessian_12_jl)
pct_diffs_34 = map(pct_diff, hessian_34_mat, hessian_34_jl)
pct_diffs_56 = map(pct_diff, hessian_56_mat, hessian_56_jl)

mean(pct_diffs_12)
mean(pct_diffs_34)
mean(pct_diffs_56)

median(pct_diffs_12)
median(pct_diffs_34)
median(pct_diffs_56)

maximum(pct_diffs_12)
maximum(pct_diffs_34)
maximum(pct_diffs_56)

pct_dict_bool = [(i, j) => pct_diffs_12[i, j] <= pct_diffs_34[i, j] <= pct_diffs_56[i, j] for i = 1:82, j = 1:82]
count(x -> x, values(pct_dict_bool))
