using Test
using HDF5
using DSGE

path = dirname(@__FILE__)

m = AnSchorfheide()
Γ0, Γ1, C, Ψ, Π = eqcond(m)
stake = 1 + 1e-6
#G1, C, impact, fmat, fwt, ywt, gev, eu, loose = gensys(Γ0, Γ1, C, Ψ, Π, stake)
G1, C, impact, eu = gensys(Γ0, Γ1, C, Ψ, Π, stake)

file = "$path/../reference/gensys.h5"
G1_exp = h5read(file, "G1_gensys")
C_exp = h5read(file,"C_gensys")
impact_exp = h5read(file, "impact")
eu_exp = h5read(file, "eu")

@testset "Check gensys outputs match reference" begin
    @test @test_matrix_approx_eq G1_exp G1
    @test @test_matrix_approx_eq C_exp C
    @test @test_matrix_approx_eq impact_exp impact

    @test isequal(eu_exp, eu)
end

G1, C, impact, eu = gensys(Γ0, Γ1, C, Ψ, Π)

# trigger indeterminacy error
m = Model1002("ss10")
Γ0, Γ1, C, Ψ, Π = eqcond(m)
Γ1[m.equilibrium_conditions[:eq_z], m.endogenous_states[:y_t]] = 1.
Γ1[m.equilibrium_conditions[:eq_z], m.endogenous_states[:y_f_t]] = 1.

@testset "Check indeterminacy zeros warning" begin
    @test_logs (:warn, "Indeterminacy: 1 loose endogenous error(s)") gensys(Γ0, Γ1, C, Ψ, Π)
end

# trigger coincident zeros
m = Model1002("ss10")
Γ0, Γ1, C, Ψ, Π = eqcond(m)
Γ1[1,:] .= 0.
Γ0[1,:] .= 0.
Ψ[1,:] .= 0.
Π[1,:] .= 0.
@testset "Check output coincident zeros warning" begin
    @test_logs (:warn, "Coincident zeros. Indeterminacy and/or nonexistence.") gensys(Γ0, Γ1, C, Ψ, Π)
    @info "After this message, a warning will be intentionally printed."
    G1, tmpC, impact, eu = gensys(Γ0, Γ1, C, Ψ, Π)
    @test eu == vec([-2, -2])
    @test isempty(G1)
    @test isempty(tmpC)
    @test isempty(impact)
end

Γ1[1,1] = 1e3
@testset "Check nonexistence warning" begin
    @test_logs (:warn, "Nonexistence: number of unstable roots exceeds number of jump variables") (:warn, "Indeterminacy: 1 loose endogenous error(s)") gensys(Γ0, Γ1, C, Ψ, Π)
end

Random.seed!(1793)
F = schur!(complex(rand(3,3)), complex(rand(3,3)))
F.T[1,1] = abs(F.S[1,1]) * (1 + 1e-3)
Γ0, Γ1, ~, ~, ~ = eqcond(m)
Fm = schur!(complex(Γ0), complex(Γ1))
@testset "Check new_div function" begin
    @test DSGE.new_div(F) == (1+5e-4)
    @test DSGE.new_div(Fm) == 1.01
end
