# Set model up
m = Model990(testing = true)
m <= Setting(:alternative_policy, DSGE.taylor99())
Γ0, Γ1, C, Ψ, Π = DSGE.taylor99_eqcond(m)

# Check policy rule
eq = m.equilibrium_conditions
endo = m.endogenous_states
@testset "Check Taylor99 rule has the right coefficients." begin
    @test Γ0[eq[:eq_mp], endo[:R_t]] == 1
    @test Γ0[eq[:eq_mp], endo[:π_a_t]] == -1.5/4
    @test Γ0[eq[:eq_mp], endo[:y_t]]  == -1/4
    @test Γ0[eq[:eq_mp], endo[:y_f_t]] == 1/4
end
