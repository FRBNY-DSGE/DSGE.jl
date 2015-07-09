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

#=
gensys_qzdiv(Γ0, Γ1, C, Ψ, Π, stake) # Runs without throwing exception
gensys_qzdiv(complex(Γ0), complex(Γ1), C, Ψ, Π, stake) # Runs without throwing exception
gensys_ordschur(Γ0, Γ1, C, Ψ, Π, stake) # Throws SingularException(1)
gensys_ordschur(complex(Γ0), complex(Γ1), C, Ψ, Π, stake) # Throws LAPACKException(1)
=#



### TEST QZ FACTORIZATION

# Matlab qz
# [AA, BB, Q, Z] = qz(G0, G1);
# alpha = diag(AA);
# beta = diag(BB);
# E = ordeig(AA, BB);
mf = MatFile("gensys/gensys_variables.mat")
AA = get_variable(mf, "AA")
BB = get_variable(mf, "BB")
Q = get_variable(mf, "Q")
Z = get_variable(mf, "Z")
alpha = get_variable(mf, "alpha")
beta = complex(get_variable(mf, "beta"))
E = get_variable(mf, "E")
close(mf)

# Julia schurfact, coercing arguments to complex
F = schurfact(complex(Γ0), complex(Γ1))
AA_schurfact, BB_schurfact, Q_schurfact, Z_schurfact = F[:S], F[:T], F[:Q]', F[:Z]


# Matlab qz vs Julia schurfact tests pass when run with `include`, not interactively
@test test_matrix_eq(AA, AA_schurfact)
@test test_matrix_eq(BB, BB_schurfact)
@test test_matrix_eq(Q, Q_schurfact)
@test test_matrix_eq(Z, Z_schurfact)
@test test_matrix_eq(alpha, F[:alpha])
@test test_matrix_eq(beta, F[:beta])
@test test_matrix_eq(E, F[:values])



### TEST QZ ORDERING

# Matlab qzdiv
# [AA_qzdiv, BB_qzdiv, Q_qzdiv, Z_qzdiv] = qzdiv(AA, BB, Q, Z);
mf = MatFile("gensys/gensys_variables.mat")
AA_qzdiv_m = get_variable(mf, "AA_qzdiv")
BB_qzdiv_m = get_variable(mf, "BB_qzdiv")
Q_qzdiv_m = get_variable(mf, "Q_qzdiv")
Z_qzdiv_m = get_variable(mf, "Z_qzdiv")
close(mf)

# Matlab ordqz
# select = abs(E) < div;
# [AA_ordqz, BB_ordqz, Q_ordqz, Z_ordqz] = ordqz(AA, BB, Q, Z, select);
mf = MatFile("gensys/gensys_variables.mat")
AA_ordqz = get_variable(mf, "AA_ordqz")
BB_ordqz = get_variable(mf, "BB_ordqz")
Q_ordqz = get_variable(mf, "Q_ordqz")
Z_ordqz = get_variable(mf, "Z_ordqz")
close(mf)

# Julia qzdiv
AA_qzdiv_j, BB_qzdiv_j, Q_qzdiv_j, Z_qzdiv_j = Gensys.qzdiv(stake, AA, BB, Q, Z)

# Julia ordschur
F2 = GeneralizedSchur(AA, BB, alpha, beta, Q', Z)
select = abs(F[:values]) .< stake
FS = ordschur(F, select)
AA_ordschur, BB_ordschur, Q_ordschur, Z_ordschur = FS[:S], FS[:T], FS[:Q]', FS[:Z]


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
@test !test_matrix_eq(AA_ordqz, AA_ordschur)
