# The given fields define the entire model structure.
# We can then concisely pass around a Model object to the remaining steps of the model
#   (solve, estimate, and forecast).
type Model990 <: AbstractModel
    spec_vars::Dict{String, Any}
    Θ::Parameters
    ind::ModelInds
    eqcond::Function
end

function Model990()
    spec = model_specifications(Model990)
    Θ = Parameters990(spec)
    ind = ModelInds(spec)
    return Model990(spec, Θ, ind, eqcond)
end

description(model::Model990) = "M990"



### Model-specific specifications
function model_specifications(T::Type{Model990})
    dict = Dict{String, Any}()
    
    # Number of anticipated policy shocks
    dict["nant"] = 6
    
    # Padding for nant
    dict["nantpad"] = 20
    
    # Number of periods back we should start incorporating zero bound expectations
    # ZLB expectations should begin in 2008 Q4
    dict["antlags"] = 24

    return dict
end





# Assign measurement equation : X_t = ZZ*S_t + DD + u_t
# where u_t = eta_t+MM* eps_t with var(eta_t) = EE
# where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
function measurement(model::Model990)
    TTT, CCC, RRR = solve(model)
    
    spec = model.spec_vars
    Θ = model.Θ
    I = model.I

    endo = I.endostates
    endo_addl = I.endostates_postgensys
    exo = I.exoshocks
    obs = I.observables

    ZZ = zeros(spec["n_observables"], spec["n_states_aug"])
    DD = zeros(spec["n_observables"], 1)
    MM = zeros(spec["n_observables"], spec["n_exoshocks"])
    EE = zeros(spec["n_observables"], spec["n_observables"])
    QQ =  zeros(spec["n_exoshocks"], spec["n_exoshocks"])


    
    ## Output growth - Quarterly!
    ZZ[obs["g_y"], endo["y_t"]] = 1.0
    ZZ[obs["g_y"], endo_addl["y_t1"]] = -1.0
    ZZ[obs["g_y"], endo["z_t"]] = 1.0
    DD[obs["g_y"]] = 100*(exp(Θ.zstar)-1)

    ## Hoursg
    ZZ[obs["hoursg"], endo["L_t"]] = 1.0
    DD[obs["hoursg"]] = Θ.Lmean.scaledvalue

    ## Labor Share/real wage growth
    ZZ[obs["g_w"], endo["w_t"]] = 1.0
    ZZ[obs["g_w"], endo_addl["w_t1"]] = -1.0
    ZZ[obs["g_w"], endo["z_t"]] = 1.0
    DD[obs["g_w"]] = 100*(exp(Θ.zstar)-1)

    ## Inflation (GDP Deflator)
    ZZ[obs["pi_gdpdef"], endo["pi_t"]] = Θ.gamm_gdpdef
    ZZ[obs["pi_gdpdef"], endo_addl["e_gdpdef"]] = 1.0
    DD[obs["pi_gdpdef"]] = 100*(Θ.pistar-1) + Θ.del_gdpdef

    ## Inflation (Core PCE)
    ZZ[obs["pi_pce"], endo["pi_t"]] = 1.0
    ZZ[obs["pi_pce"], endo_addl["e_pce"]] = 1.0
    DD[obs["pi_pce"]] = 100*(Θ.pistar-1)

    ## Nominal interest rate
    ZZ[obs["R_n"], endo["R_t"]] = 1.0
    DD[obs["R_n"]] = Θ.Rstarn

    ## Consumption Growth
    ZZ[obs["g_c"], endo["c_t"]] = 1.0
    ZZ[obs["g_c"], endo_addl["c_t1"]] = -1.0
    ZZ[obs["g_c"], endo["z_t"]] = 1.0
    DD[obs["g_c"]] = 100*(exp(Θ.zstar)-1)

    ## Investment Growth
    ZZ[obs["g_i"], endo["i_t"]] = 1.0
    ZZ[obs["g_i"], endo_addl["i_t1"]] = -1.0
    ZZ[obs["g_i"], endo["z_t"]] = 1.0
    DD[obs["g_i"]] = 100*(exp(Θ.zstar)-1)

    ## Spreads
    ZZ[obs["sprd"], endo["E_Rktil"]] = 1.0
    ZZ[obs["sprd"], endo["R_t"]] = -1.0
    DD[obs["sprd"]] = 100*log(Θ.sprd)

    ## 10 yrs infl exp
    TTT10 = (1/40)*((eye(size(TTT, 1)) - TTT)\(TTT - TTT^41))
    ZZ[obs["pi_long"], :] =  TTT10[endo["pi_t"], :]
    DD[obs["pi_long"]] = 100*(Θ.pistar-1)

    ## Long Rate
    ZZ[obs["R_long"], :] = ZZ[6, :]*TTT10
    ZZ[obs["R_long"], endo_addl["lr_t"]] = 1.0
    DD[obs["R_long"]] = Θ.Rstarn

    ## TFP
    ZZ[obs["tfp"], endo["z_t"]] = (1-Θ.alp)*Θ.modelalp_ind + 1*(1-Θ.modelalp_ind)
    ZZ[obs["tfp"], endo_addl["tfp_t"]] = 1.0
    ZZ[obs["tfp"], endo["u_t"]] = Θ.alp/( (1-Θ.alp)*(1-Θ.modelalp_ind) + 1*Θ.modelalp_ind )
    ZZ[obs["tfp"], endo_addl["u_t1"]] = -(Θ.alp/( (1-Θ.alp)*(1-Θ.modelalp_ind) + 1*Θ.modelalp_ind) )

    QQ[exo["g_sh"], exo["g_sh"]] = Θ.σ_g^2
    QQ[exo["b_sh"], exo["b_sh"]] = Θ.σ_b^2
    QQ[exo["mu_sh"], exo["mu_sh"]] = Θ.σ_mu^2
    QQ[exo["z_sh"], exo["z_sh"]] = Θ.σ_z^2
    QQ[exo["laf_sh"], exo["laf_sh"]] = Θ.σ_laf^2
    QQ[exo["law_sh"], exo["law_sh"]] = Θ.σ_law^2
    QQ[exo["rm_sh"], exo["rm_sh"]] = Θ.σ_rm^2
    QQ[exo["sigw_sh"], exo["sigw_sh"]] = Θ.σ_sigw^2
    QQ[exo["mue_sh"], exo["mue_sh"]] = Θ.σ_mue^2
    QQ[exo["gamm_sh"], exo["gamm_sh"]] = Θ.σ_gamm^2
    QQ[exo["pist_sh"], exo["pist_sh"]] = Θ.σ_pist^2
    QQ[exo["lr_sh"], exo["lr_sh"]] = Θ.σ_lr^2
    QQ[exo["zp_sh"], exo["zp_sh"]] = Θ.σ_zp^2
    QQ[exo["tfp_sh"], exo["tfp_sh"]] = Θ.σ_tfp^2
    QQ[exo["gdpdef_sh"], exo["gdpdef_sh"]] = Θ.σ_gdpdef^2
    QQ[exo["pce_sh"], exo["pce_sh"]] = Θ.σ_pce^2

    # These lines set the standard deviations for the anticipated shocks. They
    # are here no longer calibrated to the std dev of contemporaneous shocks, 
    # as we had in 904
    for i = 1:spec["nant"]
        ZZ[obs["R_n$i"], :] = ZZ[obs["R_n"], :]*(TTT^i)
        DD[obs["R_n$i"]] = Θ.Rstarn
        QQ[exo["rm_shl$i"], exo["rm_shl$i"]] = Θ.(parse("σ_rm$i"))
    end

    return ZZ, DD, QQ, EE, MM
end

