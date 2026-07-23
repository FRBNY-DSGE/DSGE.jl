sw = SmetsWouters()
Γ0, Γ1, C, Ψ, Π = eqcond(sw)
TTT_gensys, CCC_gensys, RRR_gensys, _ = gensys(Γ0, Γ1, C, Ψ, Π, 1 + 1e-6)
TTT_gensys = real(TTT_gensys)
RRR_gensys = real(RRR_gensys)
CCC_gensys = real(CCC_gensys)
TTT_aug, RRR_aug, CCC_aug = augment_states(sw, TTT_gensys, RRR_gensys, CCC_gensys)

endo     = sw.endogenous_states
endo_new = sw.endogenous_states_augmented
n_endo   = n_states(sw)           # 47
n_addl   = length(endo_new)       # 7
n_exo    = n_shocks_exogenous(sw) # 7

@testset "augment_states: output dimensions" begin
    @test size(TTT_aug) == (n_endo + n_addl, n_endo + n_addl)
    @test size(RRR_aug) == (n_endo + n_addl, n_exo)
    @test length(CCC_aug) == n_endo + n_addl
end

@testset "augment_states: original block preserved" begin
    @test TTT_aug[1:n_endo, 1:n_endo] ≈ TTT_gensys
    @test RRR_aug[1:n_endo, :] ≈ RRR_gensys
    @test CCC_aug[1:n_endo] ≈ CCC_gensys
end

@testset "augment_states: lag state assignments" begin
    @test TTT_aug[endo_new[:y_t1], endo[:y_t]] == 1.0
    @test TTT_aug[endo_new[:c_t1], endo[:c_t]] == 1.0
    @test TTT_aug[endo_new[:i_t1], endo[:i_t]] == 1.0
    @test TTT_aug[endo_new[:w_t1], endo[:w_t]] == 1.0
    @test TTT_aug[endo_new[:π_t1], endo[:π_t]] == 1.0
    @test TTT_aug[endo_new[:L_t1], endo[:L_t]]  == 1.0
end

@testset "augment_states: Et_π_t row" begin
    T2  = TTT_gensys^2
    TR  = TTT_gensys * RRR_gensys
    CTC = CCC_gensys + TTT_gensys * CCC_gensys

    @test TTT_aug[endo_new[:Et_π_t], 1:n_endo]     ≈ T2[endo[:π_t], :]
    @test TTT_aug[endo_new[:Et_π_t], n_endo+1:end] ≈ zeros(n_addl)
    @test RRR_aug[endo_new[:Et_π_t], :]            ≈ TR[endo[:π_t], :]
    @test CCC_aug[endo_new[:Et_π_t]]               ≈ CTC[endo[:π_t]]
end

@testset "augment_states: lag rows zero-initialized in RRR and CCC" begin
    lag_states = [:y_t1, :c_t1, :i_t1, :w_t1, :π_t1, :L_t1]
    for s in lag_states
        row = endo_new[s]
        @test RRR_aug[row, :] ≈ zeros(n_exo)
        @test CCC_aug[row] ≈ 0.0
    end
end

nothing
