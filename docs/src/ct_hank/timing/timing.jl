using EHANK
using BenchmarkTools

# Pre-compile
m1 = KrusellSmith()
m2 = OneAssetHANK()
~, ~, ~, ~ = solve(m1)
~, ~, ~, ~ = solve(m2)

# Time whole solution, no bypassing of checks
println("\nTiming KrusellSmith, entire solution, no bypassing checks\n")
@btime begin
    m1 = KrusellSmith()
    ~, ~, ~, ~ = solve(m1)
end

println("\nTiming OneAssetHANK, entire solution, no bypassing checks\n")
@btime begin
    m2 = OneAssetHANK()
    ~, ~, ~, ~ = solve(m2)
end

# Time whole solution, bypassing of checks
println("\nTiming KrusellSmith, entire solution but bypassing checks\n")
@btime begin
    m1 = KrusellSmith()
    test1, test2, test3, test4 = solve(m1; check_Γ0 = false)
end

println("\nTiming OneAssetHANK, entire solution but bypassing checks\n")
@btime begin
    m2 = OneAssetHANK()
    test1, test2, test3, test4 = solve(m2; check_Γ0 = false)
end

# Time whole solution, don't use sparse matrices
println("\nTiming KrusellSmith, entire solution but without sparse matrices\n")
@btime begin
    m1 = KrusellSmith()
    test1, test2, test3, test4 = solve(m1; sparse_mat = false, check_Γ0 = false)
end

println("\nTiming OneAssetHANK, entire solution but without sparse matrices\n")
@btime begin
    m2 = OneAssetHANK()
    test1, test2, test3, test4 = solve(m2; sparse_mat = false, check_Γ0 = false)
end

# Time steady state, eqcond, reduction, gensysct separately
println("\nTiming KrusellSmith, steady state only\n")
m1 = KrusellSmith()
@btime begin
    steadystate!($m1)
end

println("\nTiming OneAssetHANK, steady state only\n")
m2 = OneAssetHANK()
@btime begin
    steadystate!($m2)
end

println("\nTiming KrusellSmith, eqcond only\n")
m1 = KrusellSmith()
@btime begin
    Γ0, Γ1, Ψ, Π, C = eqcond($m1)
end

println("\nTiming OneAssetHANK, eqcond only\n")
m2 = OneAssetHANK()
@btime begin
    Γ0, Γ1, Ψ, Π, C = eqcond($m2)
end

println("\nTiming KrusellSmith, reduction only\n")
m1 = KrusellSmith()
Γ0, Γ1, Ψ, Π, C = eqcond(m1)
Γ0 = sparse(Γ0); Γ1 = sparse(Γ1); Ψ = sparse(Ψ); Π = sparse(Π); C = sparse(reshape(C, length(C), 1))
@btime begin
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_kry = krylov_reduction($m1, $Γ0, $Γ1, $Ψ, $Π, $C)
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_spl = valuef_reduction($m1, Γ0, Γ1, Ψ, Π, C)
end

println("\nTiming OneAssetHANK, reduction only\n")
m2 = OneAssetHANK()
Γ0, Γ1, Ψ, Π, C = eqcond(m2)
Γ0 = sparse(Γ0); Γ1 = sparse(Γ1); Ψ = sparse(Ψ); Π = sparse(Π); C = sparse(reshape(C, length(C), 1))
@btime begin
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_kry = krylov_reduction($m2, $Γ0, $Γ1, $Ψ, $Π, $C)
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_spl = valuef_reduction($m2, Γ0, Γ1, Ψ, Π, C)
end

println("\nTiming KrusellSmith, reduction only without sparse matrices\n")
m1 = KrusellSmith()
Γ0, Γ1, Ψ, Π, C = eqcond(m1)
@btime begin
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_kry = krylov_reduction($m1, $Γ0, $Γ1, $Ψ, $Π, $C)
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_spl = valuef_reduction($m1, Γ0, Γ1, Ψ, Π, C)
end

println("\nTiming OneAssetHANK, reduction only without sparse matrices\n")
m2 = OneAssetHANK()
Γ0, Γ1, Ψ, Π, C = eqcond(m2)
@btime begin
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_kry = krylov_reduction($m2, $Γ0, $Γ1, $Ψ, $Π, $C)
    Γ0, Γ1, Ψ, Π, C, ~, inv_basis_spl = valuef_reduction($m2, Γ0, Γ1, Ψ, Π, C)
end

println("\nTiming KrusellSmith, gensysct only\n")
m1 = KrusellSmith()
Γ0, Γ1, Ψ, Π, C = eqcond(m1)
Γ0 = sparse(Γ0); Γ1 = sparse(Γ1); Ψ = sparse(Ψ); Π = sparse(Π); C = sparse(reshape(C, length(C), 1))
Γ0, Γ1, Ψ, Π, C, ~, inv_basis_kry = krylov_reduction(m1, Γ0, Γ1, Ψ, Π, C)
Γ0, Γ1, Ψ, Π, C, ~, inv_basis_spl = valuef_reduction(m1, Γ0, Γ1, Ψ, Π, C)
Γ0 = full(Γ0); Γ1 = full(Γ1); Ψ = full(Ψ); Π = full(Π); C = full(C)
@btime begin
    TT, CC, RR, ~, ~, eu = gensysct($Γ1, $C, $Ψ, $Π, complex_decomposition = false)
end
@time begin
    TT, CC, RR, ~, ~, eu = gensysct!(Γ1, C, Ψ, Π, complex_decomposition = false)
end

println("\nTiming OneAssetHANK, gensysct only\n")
m2 = OneAssetHANK()
Γ0, Γ1, Ψ, Π, C = eqcond(m2)
Γ0 = sparse(Γ0); Γ1 = sparse(Γ1); Ψ = sparse(Ψ); Π = sparse(Π); C = sparse(reshape(C, length(C), 1))
Γ0, Γ1, Ψ, Π, C, ~, inv_basis_kry = krylov_reduction(m2, Γ0, Γ1, Ψ, Π, C)
Γ0, Γ1, Ψ, Π, C, ~, inv_basis_spl = valuef_reduction(m2, Γ0, Γ1, Ψ, Π, C)
Γ0 = full(Γ0); Γ1 = full(Γ1); Ψ = full(Ψ); Π = full(Π); C = full(C)
@btime begin
    TT, CC, RR, ~, ~, eu = gensysct($Γ1, $C, $Ψ, $Π, complex_decomposition = false)
end
@time begin
    TT, CC, RR, ~, ~, eu = gensysct!(Γ1, C, Ψ, Π, complex_decomposition = false)
end
