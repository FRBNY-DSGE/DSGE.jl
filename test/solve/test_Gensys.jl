using Base: Test, LinAlg
using MATLAB

using DSGE: Gensys
include("../util.jl")
include("../../src/solve/Gensys_versions.jl")



### TEST CALLS TO GENSYS

# Γ0, Γ1 matrices from evaluating Matlab code up to gensys call
mf = MatFile("gensys/gensys_args.mat")
Γ0 = get_variable(mf, "G0")
Γ1 = get_variable(mf, "G1")
C = get_variable(mf, "C")
Ψ = get_variable(mf, "PSI")
Π = get_variable(mf, "PIE")
stake = get_variable(mf, "div")
close(mf)

Γ0_orig, Γ1_orig, C_orig, Ψ_orig, Π_orig = copy(Γ0), copy(Γ1), copy(C), copy(Ψ), copy(Π)

#=
gensys_qzdiv(Γ0, Γ1, C, Ψ, Π, stake) # Runs without throwing exception
gensys_qzdiv(complex(Γ0), complex(Γ1), C, Ψ, Π, stake) # Runs without throwing exception
gensys_ordschur(Γ0, Γ1, C, Ψ, Π, stake) # Throws SingularException(1)
gensys_ordschur(complex(Γ0), complex(Γ1), C, Ψ, Π, stake) # Throws LAPACKException(1)
=#



### TEST QZ FACTORIZATION

# Matlab qz
mat"""
    [$AA, $BB, $Q, $Z] = qz($Γ0, $Γ1);
    $alpha = diag($AA);
    $beta = complex(diag($BB));
    $E = ordeig($AA, $BB);
"""
AA_orig, BB_orig, Q_orig, Z_orig = copy(AA), copy(BB), copy(Q), copy(Z)

# Julia schurfact, coercing arguments to complex
F = schurfact(complex(Γ0), complex(Γ1))
AA_schurfact, BB_schurfact, Q_schurfact, Z_schurfact = F[:S], F[:T], F[:Q]', F[:Z]
alpha_schurfact, beta_schurfact, E_schurfact = F[:alpha], F[:beta], F[:values]

# Matlab qz vs Julia schurfact
@test test_matrix_eq(AA, AA_schurfact)
@test test_matrix_eq(BB, BB_schurfact)
@test test_matrix_eq(Q, Q_schurfact)
@test test_matrix_eq(Z, Z_schurfact)
@test test_matrix_eq(alpha, alpha_schurfact)
@test test_matrix_eq(beta, beta_schurfact)
@test test_matrix_eq(E, E_schurfact)



### TEST QZ ORDERING

# Matlab qzdiv
# [$AA_qzdiv, $BB_qzdiv, $Q_qzdiv, $Z_qzdiv] = qzdiv($AA, $BB, $Q, $Z);
# It doesn't seem like we can call qzdiv.m using MATLAB.jl
mf = MatFile("gensys/gensys_variables.mat")
AA_qzdiv_m = get_variable(mf, "AA_qzdiv")
BB_qzdiv_m = get_variable(mf, "BB_qzdiv")
Q_qzdiv_m = get_variable(mf, "Q_qzdiv")
Z_qzdiv_m = get_variable(mf, "Z_qzdiv")
close(mf)

# Matlab ordqz
mat"""
    select = abs($E) < $stake;
   [$AA_ordqz, $BB_ordqz, $Q_ordqz, $Z_ordqz] = ordqz($AA, $BB, $Q, $Z, select);
"""

# Julia qzdiv
AA_qzdiv_j, BB_qzdiv_j, Q_qzdiv_j, Z_qzdiv_j = Gensys.qzdiv(stake, AA, BB, Q, Z)

# Julia ordschur
select = abs(F[:values]) .< stake
FS = ordschur(F, select)
AA_ordschur, BB_ordschur, Q_ordschur, Z_ordschur = FS[:S], FS[:T], FS[:Q]', FS[:Z]
alpha_ordschur, beta_ordschur, E_ordschur = FS[:alpha], FS[:beta], FS[:values]



# Matlab qzdiv and Julia qzdiv DO NOT return the same QZ ordering
@test !test_matrix_eq(AA_qzdiv_m, AA_qzdiv_j)
@test !test_matrix_eq(BB_qzdiv_m, BB_qzdiv_j)
@test !test_matrix_eq(Q_qzdiv_m, Q_qzdiv_j)
@test !test_matrix_eq(Z_qzdiv_m, Z_qzdiv_j)

# Neither do any combination of Matlab qzdiv, Matlab ordqz, Julia qzdiv, and Julia ordschur
@test !test_matrix_eq(AA_qzdiv_m, AA_ordqz)
@test !test_matrix_eq(AA_qzdiv_m, AA_ordschur)
@test !test_matrix_eq(AA_qzdiv_j, AA_ordqz)
@test !test_matrix_eq(AA_qzdiv_j, AA_ordschur)

# However, ordqz and ordschur (both call LAPACK tgsen) are much closer
# Roughly 200/4356 entries differ by > 1e-4
@test !test_matrix_eq(AA_ordqz, AA_ordschur)
@test !test_matrix_eq(BB_ordqz, BB_ordschur)
@test !test_matrix_eq(Q_ordqz, Q_ordschur)
@test test_matrix_eq(Z_ordqz, Z_ordschur)



### TEST EIGENVALUES

# alpha and beta are the diagonal elements of AA and BB after schurfact
@test test_matrix_eq(alpha_schurfact, diag(AA_schurfact))
@test test_matrix_eq(beta_schurfact, diag(BB_schurfact))
@test test_matrix_eq(E_schurfact, alpha_schurfact ./ beta_schurfact)

# No coincident zeros in alpha and beta after schurfact
coincident_zeros(x, y) = (x == 0.0 + 0.0im && y == 0.0 + 0.0im)
coincident_zeros_approx(x, y) = (abs(x) < 1e-4 && abs(y) < 1e-4)
@test !any(map(coincident_zeros, alpha_schurfact, beta_schurfact))
@test !any(map(coincident_zeros_approx, alpha_schurfact, beta_schurfact))

# Eigenvalues < stake in top left corner
@test all(x -> abs(x) < stake, E_ordschur[1:13])
@test !any(x -> abs(x) < stake, E_ordschur[14:66])

# Sort eigenvalues by abs value
E_sort = sort(E, by=abs)
E_schurfact_sort = sort(E_schurfact, by=abs)
E_ordschur_sort = sort(E_ordschur, by=abs)

# All identical in first 35 entries
@test test_matrix_eq(E_sort[1:35], E_ordschur_sort[1:35])
@test test_matrix_eq(E_schurfact_sort[1:35], E_ordschur_sort[1:35])

# All different in entries 36:51
@test all(x -> abs(x) > 1e10 && abs(x) != Inf, E_ordschur_sort[36:51])
@test all(x -> abs(x) == Inf,                  E_sort[36:66])
@test all(isnan,                               E_schurfact_sort[36:66])

# ordschur preserves decomposition
@test test_matrix_eq(Q_ordschur * Γ0 * Z_ordschur, AA_ordschur)
@test test_matrix_eq(Q_ordschur * Γ1 * Z_ordschur, BB_ordschur)



### JULIA CAN'T DIVIDE COMPLEX 1/0
@test 1.0 / 0.0 == Inf
@test isnan((1.0 + 0.0im) / (0.0 + 0.0im))



### MAKE SURE NO ARGUMENTS CHANGED DURING EVALUATION
@test Γ0 == Γ0_orig
@test Γ1 == Γ1_orig
@test C == C_orig
@test Ψ == Ψ_orig
@test Π == Π_orig

@test AA == AA_orig
@test BB == BB_orig
@test Q == Q_orig
@test Z == Z_orig
