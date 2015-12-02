using DSGE
using HDF5, Base.Test

include("../util.jl")
path = dirname(@__FILE__)

# Test hessian! in context of model
model = Model990()
model.testing = true

# Setup paths
mode = h5open(inpath(model, "user", "paramsmode.h5")) do file
    read(file, "params")
end
data = h5open(inpath(model, "data", "data_REF.h5")) do file
    read(file, "data")
end
hessian_expected = h5open(inpath(model, "user", "hessian.h5")) do file
    read(file, "hessian")
end

# Test subset of hessian elements.
para_free      = [!θ.fixed for θ in model.parameters]
para_free_inds = find(para_free)

max_free_ind = DSGE.n_hessian_test_params(model)
if max_free_ind < maximum(para_free_inds)
    para_free_inds = para_free_inds[1:max_free_ind]
end

@time hessian, _ = hessian!(model, mode, data; verbose=:none)

# The values here are liable to differ between machines, based on different
# architectures/multithreading. Hessian values tend to be large, so these larger tolerances
# are needed.
expect = hessian_expected[1:max_free_ind, 1:max_free_ind]
actual = hessian[1:max_free_ind, 1:max_free_ind]
@test_matrix_approx_eq_eps expect actual 0.1 3.0

model.testing = false

nothing
