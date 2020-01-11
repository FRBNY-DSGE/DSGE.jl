"""
```
alt_inflation()
```

Create an instance of the `AltPolicy` proposed in John Taylor's \"Discretion
versus policy rules in practice\" (1993).
"""
function alt_inflation()
    AltPolicy(:alt_inflation, alt_inflation_eqcond, alt_inflation_solve, color = RGB(1., 0., 0.))
end

function alt_inflation_eqcond(m::AbstractDSGEModel)
    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = eqcond(m)

    eq   = m.equilibrium_conditions
    endo = m.endogenous_states

    # Zero out old policy rule
    Γ0[eq[:eq_mp], :] .= 0
    Γ1[eq[:eq_mp], :] .= 0
    C[eq[:eq_mp]]     = 0
    Ψ[eq[:eq_mp], :]  .= 0
    Π[eq[:eq_mp], :]  .= 0

    # Put in new policy rule
    Γ0[eq[:eq_mp], endo[:π_t]]   = 1
    Γ0[eq[:eq_mp], endo[:π_star_t]] = -1

    return Γ0, Γ1, C, Ψ, Π
end

function alt_inflation_solve(m::AbstractDSGEModel)
    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π = alt_inflation_eqcond(m)

    # Solve model
    TTT_gensys, CCC_gensys, RRR_gensys, eu = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
    #TTT_gensys, CCC_gensys, RRR_gensys, fmat, fwt, ywt, gev, eu, loose = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
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
