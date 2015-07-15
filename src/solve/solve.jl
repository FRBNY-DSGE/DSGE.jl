using DSGE
using DSGE: AbstractModel, Gensys

# TODO: Do we need to pass in nant as an argument?
function solve(model::Model)

    Θ = model.Θ
    endo = model.I.endostates
    Γ0, Γ1, C, Ψ, Π = model.eqcond(model.spec_vars, Θ, model.I)
    T1, TC, T0, M, TZ, TY, gev, RC = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)

    TTT = real(T1)
    RRR = real(T0)
    CCC = TC

    # Some of our observables are growth rates, which is calculated as a
    # linear combination of a present and lagged state. To capture the lagged state, 
    # we assign to it an index. In addition, we also need to expand the
    # matrices of the state transition equation to accommodate the extra state.
    # In dsgesolv.m, AddTTT is appended to TTT to capture the lagged state
    # value in the current state vector, while AddRRR and AddCCC augment
    # RRR and CCC with the appropriate number of zeros.

    numAdd = 12

    # Make the Add matrices
    AddTTT = zeros(numAdd, size(TTT, 2))
    AddRRR = zeros(numAdd, size(RRR, 2))
    AddCCC = zeros(numAdd, size(CCC, 2))

    # Make the Add matrix for new states (ones that aren't just tracking old
    # lags); will be used for measurement error shocks
    AddTTT_new = zeros(numAdd, numAdd)


    
    ## AddTTT modifications

    # Track Lags
    AddTTT[1, endo["y_t"]] = 1.0
    AddTTT[2, endo["c_t"]] = 1.0
    AddTTT[3, endo["i_t"]] = 1.0
    AddTTT[4, endo["w_t"]] = 1.0
    AddTTT[5, endo["pi_t"]] = 1.0
    AddTTT[6, endo["L_t"]]  = 1.0
    AddTTT[12, endo["u_t"]] = 1.0

    # Expected inflation
    AddTTT[7, :] = (TTT^2)[endo["pi_t"], :]

    # The 8th column of AddTTT corresponds to "v_lr" which is set equal to
    # e_lr –measurements errors for the two real wage observables built in
    # as exogenous structural shocks.
    AddTTT_new[8, 8] = Θ.ρ_lr.scaledvalue
    AddTTT_new[9, 9] = Θ.ρ_tfp.scaledvalue
    AddTTT_new[10, 10] = Θ.ρ_gdpdef.scaledvalue
    AddTTT_new[11, 11] = Θ.ρ_pce.scaledvalue


    
    ## AddRRR modfications

    # Expected inflation
    AddRRR[7, :] = (TTT*RRR)[endo["pi_t"], :]

    # Measurement Error on long rate
    AddRRR[8, exo["lr_sh"]] = 1.0

    # Measurement Error on TFP
    AddRRR[9, exo["tfp_sh"]] = 1.0

    # Measurement Error on GDP Deflator
    AddRRR[10, exo["gdpdef_sh"]] = 1.0

    # Measurement Error on Core PCE
    AddRRR[11, exo["pce_sh"]] = 1.0


    
    ## AddCCC Modifications

    # Expected inflation
    AddCCC[7, :] = (CCC + TTT*CCC)[endo["pi_t"], :]

    ## Put everything together and expand TTT, RRR, CCC
    TTT = [[TTT zeros(size(TTT, 1), numAdd)]; [AddTTT AddTTT_new]]
    RRR = [RRR; AddRRR]
    CCC = [CCC; AddCCC]

    # TODO: How to deal with valid flags?
    return TTT, RRR, CCC
end
