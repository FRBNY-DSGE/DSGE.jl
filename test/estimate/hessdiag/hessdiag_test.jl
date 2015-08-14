using DSGE, MATLAB

mf = MatFile("hessdiag.mat")
hessdiag_mat = get_variable(mf, "hessdiag_mat")
hessdiag_jl = get_variable(mf, "hessdiag_jl")
close(mf)

hessian_12_mat = zeros(82, 82)
hessian_12_jl = zeros(82, 82)
hessian_34_mat = zeros(82, 82)
hessian_34_jl = zeros(82, 82)
for i = 1:82
    for j = 1:82
        hessian_12_mat[i, j] = -(hessdiag_mat[i, j, 1] + hessdiag_mat[i, j, 2])/2
        hessian_12_jl[i, j] = -(hessdiag_jl[i, j, 1] + hessdiag_jl[i, j, 2])/2
        hessian_34_mat[i, j] = -(hessdiag_mat[i, j, 3] + hessdiag_mat[i, j, 4])/2
        hessian_34_jl[i, j] = -(hessdiag_jl[i, j, 3] + hessdiag_jl[i, j, 4])/2
    end
end

function pct_diff(x, y)
    abs_diff = abs(abs(x) - abs(y))
    if x != 0
        return abs_diff/abs(x)
    elseif y != 0
        return abs_diff/abs(y)
    else
        return 0
    end
end

pct_diffs_12 = map(pct_diff, hessian_12_mat, hessian_12_jl)
pct_diffs_34 = map(pct_diff, hessian_34_mat, hessian_34_jl)
