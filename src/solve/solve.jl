using ..AbstractModel, ..Gensys

# TODO: Do we need to pass in nant as an argument?
function solve(model::Model)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π = model.eqcond(model.spec_vars, model.Θ, model.I)

    # Solve model
    TTT_gensys, CCC_gensys, RRR_gensys = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
    TTT_gensys = real(TTT_gensys)
    RRR_gensys = real(RRR_gensys)

    # Augment states
    TTT, CCC, RRR = augment_states(model, TTT_gensys, CCC_gensys, RRR_gensys)

    return TTT, CCC, RRR
end



# Some of our observables are growth rates, which is calculated as a
# linear combination of a present and lagged state. To capture the lagged state, 
# we assign to it an index. In addition, we also need to expand the
# matrices of the state transition equation to accommodate the extra state.
# In dsgesolv.m, AddTTT is appended to TTT to capture the lagged state
# value in the current state vector, while AddRRR and AddCCC augment
# RRR and CCC with the appropriate number of zeros.

# These additional states are added after the model is solved to reduce the load on gensys
function augment_states{T<:FloatingPoint}(model::Model, TTT::Array{T, 2}, CCC::Array{T, 2}, RRR::Array{T, 2})
    Θ = model.Θ
    endo = model.I.endostates
    endo_addl = model.I.endostates_postgensys

    n_endo = length(endo)
    n_exo = length(exo)
    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo) == size(RRR)
    @assert (n_endo, 1) == size(CCC)

    # Initialize augmented matrices
    numAdd = 12
    TTT_aug = zeros(n_endo + numAdd, n_endo + numAdd)
    TTT_aug[1:n_endo, 1:n_endo] = TTT
    RRR_aug = [RRR; zeros(numAdd, n_exo)]
    CCC_aug = [CCC; zeros(numAdd, 1)]


    
    ### TTT modifications

    # Track Lags
    TTT_aug[endo_addl["y_t1"], endo["y_t"]] = 1.0
    TTT_aug[endo_addl["c_t1"], endo["c_t"]] = 1.0
    TTT_aug[endo_addl["i_t1"], endo["i_t"]] = 1.0
    TTT_aug[endo_addl["w_t1"], endo["w_t"]] = 1.0
    TTT_aug[endo_addl["pi_t1"], endo["pi_t"]] = 1.0
    TTT_aug[endo_addl["w_t1"], endo["L_t"]]  = 1.0
    TTT_aug[endo_addl["u_t1"], endo["u_t"]] = 1.0

    # Expected inflation
    TTT_aug[endo_addl["Et_pi_t"], 1:n_endo] = (TTT^2)[endo["pi_t"], :]

    # The 8th column of AddTTT corresponds to "v_lr" which is set equal to
    # e_lr –measurements errors for the two real wage observables built in
    # as exogenous structural shocks.
    TTT_aug[endo_addl["lr_t"], endo_addl["lr_t"]] = Θ.ρ_lr.scaledvalue
    TTT_aug[endo_addl["tfp_t"], endo_addl["tfp_t"]] = Θ.ρ_tfp.scaledvalue
    TTT_aug[endo_addl["e_gdpdef"], endo_addl["e_gdpdef"]] = Θ.ρ_gdpdef.scaledvalue
    TTT_aug[endo_addl["e_gdpdef"], endo_addl["e_pce"]] = Θ.ρ_pce.scaledvalue


    
    ### CCC Modifications

    # Expected inflation
    CCC_aug[endo_addl["Et_pi_t"], :] = (CCC + TTT*CCC)[endo["pi_t"], :]


    
    ### RRR modfications

    # Expected inflation
    RRR_aug[endo_addl["Et_pi_t"], :] = (TTT*RRR)[endo["pi_t"], :]

    # Measurement Error on long rate
    RRR_aug[endo_addl["lr_t"], exo["lr_sh"]] = 1.0

    # Measurement Error on TFP
    RRR_aug[endo_addl["tfp_t"], exo["tfp_sh"]] = 1.0

    # Measurement Error on GDP Deflator
    RRR_aug[endo_addl["e_gdpdef"], exo["gdpdef_sh"]] = 1.0

    # Measurement Error on Core PCE
    RRR_aug[endo_addl["e_pce"], exo["pce_sh"]] = 1.0



    # TODO: How to deal with valid flags?
    return TTT_aug, CCC_aug, RRR_aug
end
