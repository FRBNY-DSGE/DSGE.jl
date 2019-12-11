using DSGE, SparseArrays, MAT, FileIO

# Script for doing all that MATLAB does
#m = TwoAssetHANK()
#steadystate!(m)
reduceDistribution = true #get_setting(m, :reduceDistribution)
reduceDist_hor     = 100 #get_setting(m, :reduceDist_hor)
reduceV            = 1 #get_setting(m, :reduceV)

nVars = 127504 #get_setting(m, :nVars)
n_v = 60000 #get_setting(m, :n_v)
n_g = 67499 #get_setting(m, :n_g)
n_p = 4 #get_setting(m, :n_p)
n_Z = 1 #get_setting(m, :n_Z)

#Γ0, Γ1, Ψ, Π C = eqcond(m)
test_out = load("data/eqcond_output_matlab.jld2")
Γ0 = test_out["g0"]
Γ1 = test_out["g1"]
Ψ = test_out["psi"]
Π = test_out["pi"]
C = spzeros(Float64, nVars)

# Distribution Reduction
base_mat = my_speye(nVars)
inv_base = my_speye(nVars)
state_red = my_speye(nVars-n_p)
inv_state_red = my_speye(nVars-n_p)

println("=> Reducing Value Function Distribution")
@show reduceDistribution
#error()
@time if reduceDistribution

   state_red, inv_state_red, n_g = stateSpaceReduction(Γ0, Γ1, n_v, n_g, n_p, n_Z, reduceDist_hor)

else
    Γ0, Γ1, C, Π, Ψ, base_mat, inv_base = solve_static_conditions(Γ0, Γ1, Π, Ψ, C)
    # Confirming these do the same thing
    base_mat2, inv_base2, g1, c, pi, psi  = cleanG0Sparse(Γ0, Γ1, C, Π, Ψ)
    @show base_mat2 ≈ base_mat, inv_base ≈ inv_base2, g1 ≈ Γ1, Π ≈ pi, Ψ ≈ psi
end

# Value Function Reduction
from_spline = my_speye(n_g + n_v + n_Z)
to_spline   = my_speye(n_g + n_v + n_Z)
n_splined   = n_v
n_rest      = n_g + n_Z

if reduceV
    #to_spline, from_spline, n_splined, n_rest = reduce
end
