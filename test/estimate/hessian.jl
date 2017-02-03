using DSGE
using HDF5, Base.Test

include("../util.jl")
path = dirname(@__FILE__)

# # Test hessian! in context of model
m = AnSchorfheide()
m.testing = true

# Setup paths

data = h5read("$path/../reference/hessian.h5","data")
mode = h5read("$path/../reference/hessian.h5","paramsmode")
hessian_expected = h5read("$path/../reference/hessian.h5","hessian")

# Test subset of hessian elements.
para_free      = [!θ.fixed for θ in m.parameters]
para_free_inds = find(para_free)

max_free_ind = DSGE.n_hessian_test_params(m)
if max_free_ind < maximum(para_free_inds)
    para_free_inds = para_free_inds[1:max_free_ind]
end

@time hessian, _ = hessian!(m, mode, data; verbose=:none)

# The values here are liable to differ between machines, based on different
# architectures/multithreading. Hessian values tend to be large, so these larger tolerances
# are needed.
expect = hessian_expected[1:max_free_ind, 1:max_free_ind]
actual = hessian[1:max_free_ind, 1:max_free_ind]
@test_matrix_approx_eq_eps expect actual 0.1 3.0

m.testing = false

nothing
