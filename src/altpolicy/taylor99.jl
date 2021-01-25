"""
```
taylor99()
```

Create an instance of the `AltPolicy` proposed in John Taylor's \"A Historical
Analysis of Monetary Policy Rules\" (1999).
"""
function taylor99()
    AltPolicy(:taylor99, taylor99_eqcond, taylor99_solve, color = RGB(0., 0., 1.))
end

function taylor99_eqcond(m::AbstractDSGEModel, reg::Int = 1)
    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = eqcond(m, reg)
    eq   = m.equilibrium_conditions
    endo = m.endogenous_states

    # Zero out old policy rule
    Γ0[eq[:eq_mp], :] .= 0
    Γ1[eq[:eq_mp], :] .= 0
    C[eq[:eq_mp]]     = 0
    Ψ[eq[:eq_mp], :]  .= 0
    Π[eq[:eq_mp], :]  .= 0

    # Put in new policy rule
    Γ0[eq[:eq_mp], endo[:R_t]]   = 1
    Γ0[eq[:eq_mp], endo[:π_a_t]] = -1.5/4
    Γ0[eq[:eq_mp], endo[:y_t]]   = -1/4
    Γ0[eq[:eq_mp], endo[:y_f_t]] = 1/4

    return Γ0, Γ1, C, Ψ, Π
end

function taylor99_solve(m::AbstractDSGEModel)
    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = taylor99_eqcond(m)

    # Solve model
    #TTT_gensys, CCC_gensys, RRR_gensys, fmat, fwt, ywt, gev, eu, loose = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
    TTT_gensys, CCC_gensys, RRR_gensys, eu  = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
    if !((eu[1] == 1) & (eu[2] == 1))
        throw(GensysError("Gensys does not give existence"))
    end
    TTT_gensys = real(TTT_gensys)
    RRR_gensys = real(RRR_gensys)
    CCC_gensys = reshape(CCC_gensys, size(CCC_gensys, 1))

    # Augment states
    TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)

    return TTT, RRR, CCC
end
