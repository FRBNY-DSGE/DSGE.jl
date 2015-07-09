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
gensys_ordschur(Γ0, Γ1, C, Ψ, Π, stake) # Throws a SingularException
=#



### TEST QZ FACTORIZATION

# Matlab qz
# [AA, BB, Q, Z] = qz(G0, G1);
# alpha = diag(AA);
# beta = diag(BB);
mf = MatFile("gensys/gensys_variables.mat")
AA = get_variable(mf, "AA")
BB = get_variable(mf, "BB")
Q = get_variable(mf, "Q")
Z = get_variable(mf, "Z")
alpha = get_variable(mf, "alpha")
beta = complex(get_variable(mf, "beta"))
close(mf)

# Julia schurfact
F = schurfact(complex(Γ0), complex(Γ1))
AA_schurfact, BB_schurfact, Q_schurfact, Z_schurfact = F[:S], F[:T], F[:Q]', F[:Z]

# Matlab qz vs Julia schurfact
# These all pass.
@test test_matrix_eq(AA, AA_schurfact)
@test test_matrix_eq(BB, BB_schurfact)
@test test_matrix_eq(Q, Q_schurfact)
@test test_matrix_eq(Z, Z_schurfact)

### TEST QZ ORDERING

# Matlab qzdiv
# [AA_qzdiv, BB_qzdiv, Q_qzdiv, Z_qzdiv] = qzdiv(AA, BB, Q, Z);
mf = MatFile("gensys/gensys_variables.mat")
AA_qzdiv = get_variable(mf, "AA_qzdiv")
BB_qzdiv = get_variable(mf, "BB_qzdiv")
Q_qzdiv = get_variable(mf, "Q_qzdiv")
Z_qzdiv = get_variable(mf, "Z_qzdiv")
close(mf)

# Matlab ordqz
# E = ordeig(AA, BB);
# select = abs(E) < div;
# [AA_ordqz, BB_ordqz, Q_ordqz, Z_ordqz] = ordqz(AA, BB, Q, Z, select);
mf = MatFile("gensys/gensys_variables.mat")
AA_ordqz = get_variable(mf, "AA_ordqz")
BB_ordqz = get_variable(mf, "BB_ordqz")
Q_ordqz = get_variable(mf, "Q_ordqz")
Z_ordqz = get_variable(mf, "Z_ordqz")
close(mf)

# Julia qzdiv
AA_qzdiv_j, BB_qzdiv_j, Q_qzdiv_j, Z_qzdiv_j = Gensys.qzdiv(stake, copy(AA), copy(BB), copy(Q), copy(Z))
AA_qzdiv_j_CC, BB_qzdiv_j_CC, Q_qzdiv_j_CC, Z_qzdiv_j_CC = qzdiv_CC(stake, copy(AA), copy(BB), copy(Q), copy(Z))

# Julia ordschur
F = GeneralizedSchur(copy(AA), copy(BB), alpha, beta, copy(Q)', copy(Z))
select = abs(F[:values]) .< stake
FS = ordschur(F, select)
AA_ordschur, BB_ordschur, Q_ordschur, Z_ordschur = FS[:S], FS[:T], FS[:Q]', FS[:Z]

#=
# Matlab qzdiv vs Julia qzdiv
# These don't pass
@test test_matrix_eq(AA_qzdiv, AA_qzdiv_j)
@test test_matrix_eq(BB_qzdiv, BB_qzdiv_j)
@test test_matrix_eq(Q_qzdiv, Q_qzdiv_j)
@test test_matrix_eq(Z_qzdiv, Z_qzdiv_j)
=#

#=
# Matlab qzdiv vs Julia qzdiv_CC
# I could have sworn these passed at one point, but now they don't??
@test test_matrix_eq(AA_qzdiv, AA_qzdiv_j_CC)
@test test_matrix_eq(BB_qzdiv, BB_qzdiv_j_CC)
@test test_matrix_eq(Q_qzdiv, Q_qzdiv_j_CC)
@test test_matrix_eq(Z_qzdiv, Z_qzdiv_j_CC)
=#

#=
# None of these pass
@test test_matrix_eq(AA_qzdiv, AA_ordschur)
@test test_matrix_eq(AA_ordqz, AA_ordschur)
@test test_matrix_eq(AA_qzdiv, AA_ordqz)
=#
