using EHANK, BenchmarkTools

# Instantiate model objects, prep for differentiation
m1 = KrusellSmith()
m2 = OneAssetHANK()

x1 = zeros(Float64, 2 * n_states(m1) + n_shocks_expectational(m1) + n_shocks_exogenous(m1))
x2 = zeros(Float64, 2 * n_states(m2) + n_shocks_expectational(m2) + n_shocks_exogenous(m2))

function f1{T<:Real}(x::Vector{T});
    v_residual = EHANK.get_residuals(x, m1)
    return vcat(values(v_residual)...)
end

function f2{T<:Real}(x::Vector{T});
    v_residual = EHANK.get_residuals(x, m2)
    return vcat(values(v_residual)...)
end

# Time differentiation step
derivs1 = zeros(Float64, n_states(m1), length(x1))
derivs2 = zeros(Float64, n_states(m2), length(x2))
eqcond(m1); # make sure compiled
eqcond(m2);

println("Timing KrusellSmith, get_residuals")
@btime begin
    out = f1(x1);
end

println("Timing OneAssetHANK, get_residuals")
@btime begin
    out = f2(x2);
end

println("Timing KrusellSmith, jacobian")
@btime begin
    derivs1[:,:] = ForwardDiff.jacobian(f1,x1)::Matrix{Float64}
end

println("Timing OneAssetHANK, jacobian")
@btime begin
    derivs2[:,:] = ForwardDiff.jacobian(f2,x2)::Matrix{Float64}
end
