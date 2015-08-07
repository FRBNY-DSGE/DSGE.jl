# The given fields define the entire model structure.
# We can then concisely pass around a Model object to the remaining steps of the model
#   (solve, estimate, and forecast).
type Model990 <: AbstractModel
    spec::Dict{String, Any}
    Θ::Parameters
    ind::ModelInds
end

function Model990()
    spec = model_specifications(Model990)
    Θ = Parameters990(spec)
    ind = ModelInds(spec)
    return Model990(spec, Θ, ind)
end





# Model-specific specifications
function model_specifications(T::Type{Model990})
    spec = Dict{String, Any}()

    
    
    ### Model initialization
    
    # Number of anticipated policy shocks
    spec["nant"] = 6

    # Padding for nant
    spec["nantpad"] = 20

    # Number of periods back we should start incorporating zero bound expectations
    # ZLB expectations should begin in 2008 Q4
    spec["antlags"] = 24

    # TODO: This should be set when the data are read in
    # Number of presample periods
    spec["n_presample_periods"] = 2



    ### Estimation

    # Whether or not to call csminwel
    spec["reoptimize"] = false

    # Whether or not to recalculate hessian
    spec["recalculate_hessian"] = false
    
    # Number of draws per posterior simulation block
    spec["n_sim"] = 10000

    # Number of blocks
    spec["n_blocks"] = 22

    # Save every ntimes-th draw that is accepted
    spec["n_times"] = 5

    
    
    return spec
end





# Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
# Using the assigned states and equations in modelinds.jl, coefficients are specified in the
#   proper positions.

# Γ0 (n states x n states) holds coefficients of current time states.
# Γ1 (n states x n states) holds coefficients of lagged states.
# C (n states x 1) is a vector of constants
# Ψ (n states x n exogenous shocks) holds coefficients of iid shocks.
# Π (n states x n expectational states) holds coefficients of expectational states.
function eqcond(model::AbstractModel)
    spec = model.spec
    Θ = model.Θ
    ind = model.ind

    endo = ind.endostates
    exo = ind.exoshocks
    ex = ind.expshocks
    eq = ind.eqconds

    Γ0 = zeros(spec["n_states"], spec["n_states"])
    Γ1 = zeros(spec["n_states"], spec["n_states"])
    C =  zeros(spec["n_states"], 1)
    Ψ = zeros(spec["n_states"], spec["n_exoshocks"])
    Π =  zeros(spec["n_states"], spec["n_expshocks"])



    ### ENDOGENOUS STATES ###

    ### 1. Consumption Euler Equation

    # Sticky prices and wages
    Γ0[eq["euler"], endo["c_t"]] = 1.
    Γ0[eq["euler"], endo["R_t"]] = (1 - Θ.h*exp(-Θ.zstar))/(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))
    Γ0[eq["euler"], endo["b_t"]] = -1.
    Γ0[eq["euler"], endo["E_pi"]] = -(1 - Θ.h*exp(-Θ.zstar))/(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))
    Γ0[eq["euler"], endo["z_t"]] = (Θ.h*exp(-Θ.zstar))/(1 + Θ.h*exp(-Θ.zstar))
    Γ0[eq["euler"], endo["E_c"]] = -1/(1 + Θ.h*exp(-Θ.zstar))
    Γ0[eq["euler"], endo["E_z"]] = -1/(1 + Θ.h*exp(-Θ.zstar))
    Γ0[eq["euler"], endo["L_t"]] = -(Θ.sigmac - 1)*Θ.wl_c/(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))
    Γ0[eq["euler"], endo["E_L"]] = (Θ.sigmac - 1)*Θ.wl_c/(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))
    Γ1[eq["euler"], endo["c_t"]] = (Θ.h*exp(-Θ.zstar))/(1 + Θ.h*exp(-Θ.zstar))

    # Flexible prices and wages
    Γ0[eq["euler_f"], endo["c_f_t"]] = 1.
    Γ0[eq["euler_f"], endo["r_f_t"]] = (1 - Θ.h*exp(-Θ.zstar))/(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))
    Γ0[eq["euler_f"], endo["b_t"]] = -1.
    Γ0[eq["euler_f"], endo["z_t"]] =   (Θ.h*exp(-Θ.zstar))/(1 + Θ.h*exp(-Θ.zstar))
    Γ0[eq["euler_f"], endo["E_c_f"]] = -1/(1 + Θ.h*exp(-Θ.zstar))
    Γ0[eq["euler_f"], endo["E_z"]] = -1/(1 + Θ.h*exp(-Θ.zstar))
    Γ0[eq["euler_f"], endo["L_f_t"]] = -(Θ.sigmac - 1)*Θ.wl_c/(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))
    Γ0[eq["euler_f"], endo["E_L_f"]] = (Θ.sigmac - 1)*Θ.wl_c/(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))
    Γ1[eq["euler_f"], endo["c_f_t"]] = (Θ.h*exp(-Θ.zstar))/(1 + Θ.h*exp(-Θ.zstar))



    ### 2. Investment Euler Equation

    # Sticky prices and wages
    Γ0[eq["inv"], endo["qk_t"]] = -1/(Θ.s2*exp(2.*Θ.zstar)*(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)))
    Γ0[eq["inv"], endo["i_t"]] = 1.
    Γ0[eq["inv"], endo["z_t"]] = 1/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ1[eq["inv"], endo["i_t"]] = 1/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["inv"], endo["E_i"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["inv"], endo["E_z"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["inv"], endo["mu_t"]] = -1.

    # Flexible prices and wages
    Γ0[eq["inv_f"], endo["qk_f_t"]] = -1/(Θ.s2*exp(2*Θ.zstar)*(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)))
    Γ0[eq["inv_f"], endo["i_f_t"]] = 1.
    Γ0[eq["inv_f"], endo["z_t"]] = 1/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ1[eq["inv_f"], endo["i_f_t"]] = 1/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["inv_f"], endo["E_i_f"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["inv_f"], endo["E_z"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["inv_f"], endo["mu_t"]] = -1.



    ### 3. Financial Friction Block

    # Return to capital
    # Sticky prices and wages
    Γ0[eq["capval"], endo["Rktil_t"]] = 1.
    Γ0[eq["capval"], endo["pi_t"]] = -1.
    Γ0[eq["capval"], endo["rk_t"]] = -Θ.rkstar/(1 + Θ.rkstar - Θ.del)
    Γ0[eq["capval"], endo["qk_t"]] = -(1 - Θ.del)/(1 + Θ.rkstar - Θ.del)
    Γ1[eq["capval"], endo["qk_t"]] = -1.

    # Spreads
    # Sticky prices and wages
    Γ0[eq["spread"], endo["E_Rktil"]] = 1.
    Γ0[eq["spread"], endo["R_t"]] = -1.
    Γ0[eq["spread"], endo["b_t"]] = (Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))/(1 - Θ.h*exp(-Θ.zstar))
    Γ0[eq["spread"], endo["qk_t"]] = -Θ.zeta_spb
    Γ0[eq["spread"], endo["kbar_t"]] = -Θ.zeta_spb
    Γ0[eq["spread"], endo["n_t"]] = Θ.zeta_spb
    Γ0[eq["spread"], endo["sigw_t"]] = -1.
    Γ0[eq["spread"], endo["mue_t"]] = -1.

    # n evol
    # Sticky prices and wages
    Γ0[eq["nevol"], endo["n_t"]] = 1.
    Γ0[eq["nevol"], endo["gamm_t"]] = -1.
    Γ0[eq["nevol"], endo["z_t"]] = Θ.gammstar*Θ.vstar/Θ.nstar
    Γ0[eq["nevol"], endo["Rktil_t"]] = -Θ.zeta_nRk
    Γ0[eq["nevol"], endo["pi_t"]] = (Θ.zeta_nRk - Θ.zeta_nR)
    Γ1[eq["nevol"], endo["sigw_t"]] = -Θ.zeta_nsigw/Θ.zeta_spsigw
    Γ1[eq["nevol"], endo["mue_t"]] = -Θ.zeta_nmue/Θ.zeta_spmue
    Γ1[eq["nevol"], endo["qk_t"]] = Θ.zeta_nqk
    Γ1[eq["nevol"], endo["kbar_t"]] = Θ.zeta_nqk
    Γ1[eq["nevol"], endo["n_t"]] = Θ.zeta_nn
    Γ1[eq["nevol"], endo["R_t"]] = -Θ.zeta_nR
    Γ1[eq["nevol"], endo["b_t"]] = -Θ.zeta_nR

    # Flexible prices and wages - ASSUME NO FINANCIAL FRICTIONS
    Γ0[eq["capval_f"], endo["E_rk_f"]] = -Θ.rkstar/(1 + Θ.rkstar - Θ.del)
    Γ0[eq["capval_f"], endo["E_qk_f"]] = -(1 - Θ.del)/(1 + Θ.rkstar - Θ.del)
    Γ0[eq["capval_f"], endo["qk_f_t"]] = 1.
    Γ0[eq["capval_f"], endo["r_f_t"]] = 1.
    Γ0[eq["capval_f"], endo["b_t"]] = -(Θ.sigmac*(1 + Θ.h*exp(-Θ.zstar)))/(1 - Θ.h*exp(-Θ.zstar))



    ### 4. Aggregate Production Function

    # Sticky prices and wages
    Γ0[eq["output"], endo["y_t"]] =  1.
    Γ0[eq["output"], endo["k_t"]] = -Θ.Bigphi*Θ.alp
    Γ0[eq["output"], endo["L_t"]] = -Θ.Bigphi*(1 - Θ.alp)

    # Flexible prices and wages
    Γ0[eq["output_f"], endo["y_f_t"]] =  1.
    Γ0[eq["output_f"], endo["k_f_t"]] = -Θ.Bigphi*Θ.alp
    Γ0[eq["output_f"], endo["L_f_t"]] = -Θ.Bigphi*(1 - Θ.alp)



    ### 5. Capital Utilization

    # Sticky prices and wages
    Γ0[eq["caputl"], endo["k_t"]] =  1.
    Γ1[eq["caputl"], endo["kbar_t"]] =  1.
    Γ0[eq["caputl"], endo["z_t"]] = 1.
    Γ0[eq["caputl"], endo["u_t"]] = -1.

    # Flexible prices and wages
    Γ0[eq["caputl_f"], endo["k_f_t"]] =  1.
    Γ1[eq["caputl_f"], endo["kbar_f_t"]] =  1.
    Γ0[eq["caputl_f"], endo["z_t"]] = 1.
    Γ0[eq["caputl_f"], endo["u_f_t"]] = -1.



    ### 6. Rental Rate of Capital

    # Sticky prices and wages
    Γ0[eq["capsrv"], endo["u_t"]] = 1.
    Γ0[eq["capsrv"], endo["rk_t"]] = -(1 - Θ.ppsi)/Θ.ppsi

    # Flexible prices and wages
    Γ0[eq["capsrv_f"], endo["u_f_t"]] = 1.
    Γ0[eq["capsrv_f"], endo["rk_f_t"]] = -(1 - Θ.ppsi)/Θ.ppsi



    ### 7. Evolution of Capital

    # Sticky prices and wages
    Γ0[eq["capev"], endo["kbar_t"]] = 1.
    Γ1[eq["capev"], endo["kbar_t"]] = 1 - Θ.istar/Θ.kbarstar
    Γ0[eq["capev"], endo["z_t"]] = 1 - Θ.istar/Θ.kbarstar
    Γ0[eq["capev"], endo["i_t"]] = -Θ.istar/Θ.kbarstar
    Γ0[eq["capev"], endo["mu_t"]] = -Θ.istar*Θ.s2*exp(2*Θ.zstar)*(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))/Θ.kbarstar

    # Flexible prices and wages
    Γ0[eq["capev_f"], endo["kbar_f_t"]] = 1.
    Γ1[eq["capev_f"], endo["kbar_f_t"]] = 1 - Θ.istar/Θ.kbarstar
    Γ0[eq["capev_f"], endo["z_t"]] = 1 - Θ.istar/Θ.kbarstar
    Γ0[eq["capev_f"], endo["i_f_t"]] = -Θ.istar/Θ.kbarstar
    Γ0[eq["capev_f"], endo["mu_t"]] = -Θ.istar*Θ.s2*exp(2*Θ.zstar)*(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))/Θ.kbarstar



    ### 8. Price Markup

    # Sticky prices and wages
    Γ0[eq["mkupp"], endo["mc_t"]] = 1
    Γ0[eq["mkupp"], endo["w_t"]] =  -1.
    Γ0[eq["mkupp"], endo["L_t"]] =  -Θ.alp
    Γ0[eq["mkupp"], endo["k_t"]] =  Θ.alp

    # Flexible prices and wages
    Γ0[eq["mkupp_f"], endo["w_f_t"]] = 1.
    Γ0[eq["mkupp_f"], endo["L_f_t"]] =  Θ.alp
    Γ0[eq["mkupp_f"], endo["k_f_t"]] =  -Θ.alp



    ### 9. Phillips Curve

    # Sticky prices and wages
    Γ0[eq["phlps"], endo["pi_t"]] = 1.
    Γ0[eq["phlps"], endo["mc_t"]] =  -((1 - Θ.zeta_p*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))*(1 - Θ.zeta_p))/(Θ.zeta_p*((Θ.Bigphi- 1)*Θ.epsp + 1))/(1 + Θ.iota_p*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ1[eq["phlps"], endo["pi_t"]] = Θ.iota_p/(1 + Θ.iota_p*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["phlps"], endo["E_pi"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.iota_p*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    # Comment out for counterfactual with no price mark up shock
    Γ0[eq["phlps"], endo["laf_t"]] = -(1 + Θ.iota_p*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))/(1 + Θ.iota_p*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))

    # Flexible prices and wages not necessary



    ### 10. Rental Rate of Capital

    # Sticky prices and wages
    Γ0[eq["caprnt"], endo["rk_t"]] = 1.
    Γ0[eq["caprnt"], endo["k_t"]] = 1.
    Γ0[eq["caprnt"], endo["L_t"]] = -1.
    Γ0[eq["caprnt"], endo["w_t"]] = -1.

    # Flexible prices and wages
    Γ0[eq["caprnt_f"], endo["rk_f_t"]] = 1.
    Γ0[eq["caprnt_f"], endo["k_f_t"]] = 1.
    Γ0[eq["caprnt_f"], endo["L_f_t"]] = -1.
    Γ0[eq["caprnt_f"], endo["w_f_t"]] = -1.



    ### 11. Marginal Substitution

    # Sticky prices and wages
    Γ0[eq["msub"], endo["muw_t"]] = 1.
    Γ0[eq["msub"], endo["L_t"]] = Θ.nu_l
    Γ0[eq["msub"], endo["c_t"]] = 1/(1 - Θ.h*exp(-Θ.zstar))
    Γ1[eq["msub"], endo["c_t"]] = Θ.h*exp(-Θ.zstar)/(1 - Θ.h*exp(-Θ.zstar))
    Γ0[eq["msub"], endo["z_t"]] = Θ.h*exp(-Θ.zstar) /(1 - Θ.h*exp(-Θ.zstar))
    Γ0[eq["msub"], endo["w_t"]] = -1.

    # Flexible prices and wages
    Γ0[eq["msub_f"], endo["w_f_t"]] = -1.
    Γ0[eq["msub_f"], endo["L_f_t"]] = Θ.nu_l
    Γ0[eq["msub_f"], endo["c_f_t"]] = 1/(1 - Θ.h*exp(-Θ.zstar))
    Γ1[eq["msub_f"], endo["c_f_t"]] = Θ.h*exp(-Θ.zstar)/(1 - Θ.h*exp(-Θ.zstar))
    Γ0[eq["msub_f"], endo["z_t"]] = Θ.h*exp(-Θ.zstar)/(1 - Θ.h*exp(-Θ.zstar))



    ### 12. Evolution of Wages

    # Sticky prices and wages
    Γ0[eq["wage"], endo["w_t"]] = 1
    Γ0[eq["wage"], endo["muw_t"]] = (1 - Θ.zeta_w*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))*(1 - Θ.zeta_w)/(Θ.zeta_w*((Θ.law - 1)*Θ.epsw + 1))/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["wage"], endo["pi_t"]] = (1 + Θ.iota_w*Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ1[eq["wage"], endo["w_t"]] = 1/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["wage"], endo["z_t"]] = 1/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ1[eq["wage"], endo["pi_t"]] = Θ.iota_w/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["wage"], endo["E_w"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["wage"], endo["E_z"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["wage"], endo["E_pi"]] = -Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar)/(1 + Θ.bet*exp((1 - Θ.sigmac)*Θ.zstar))
    Γ0[eq["wage"], endo["law_t"]] = -1.

    # Flexible prices and wages not necessary



    ### 13. Monetary Policy Rule

    # Sticky prices and wages
    Γ0[eq["mp"], endo["R_t"]] = 1.
    Γ1[eq["mp"], endo["R_t"]] = Θ.rho
    Γ0[eq["mp"], endo["pi_t"]] = -(1 - Θ.rho)*Θ.psi1
    Γ0[eq["mp"], endo["pist_t"]] = (1 - Θ.rho)*Θ.psi1
    Γ0[eq["mp"], endo["y_t"]] = -(1 - Θ.rho)*Θ.psi2 - Θ.psi3
    Γ0[eq["mp"], endo["y_f_t"]] = (1 - Θ.rho)*Θ.psi2 + Θ.psi3
    Γ1[eq["mp"], endo["y_t"]] = -Θ.psi3
    Γ1[eq["mp"], endo["y_f_t"]] = Θ.psi3
    Γ0[eq["mp"], endo["rm_t"]] = -1.

    # Flexible prices and wages not necessary



    ### 14. Resource Constraint

    # Sticky prices and wages
    Γ0[eq["res"], endo["y_t"]] = 1.
    Γ0[eq["res"], endo["g_t"]] = -Θ.gstar
    Γ0[eq["res"], endo["c_t"]] = -Θ.cstar/Θ.ystar
    Γ0[eq["res"], endo["i_t"]] = -Θ.istar/Θ.ystar
    Γ0[eq["res"], endo["u_t"]] = -Θ.rkstar*Θ.kstar/Θ.ystar

    # Flexible prices and wages
    Γ0[eq["res_f"], endo["y_f_t"]] = 1.
    Γ0[eq["res_f"], endo["g_t"]] = -Θ.gstar
    Γ0[eq["res_f"], endo["c_f_t"]] = -Θ.cstar/Θ.ystar
    Γ0[eq["res_f"], endo["i_f_t"]] = -Θ.istar/Θ.ystar
    Γ0[eq["res_f"], endo["u_f_t"]] = -Θ.rkstar*Θ.kstar/Θ.ystar



    ### 15. Extra States
    # These aren't strictly necessary, but they track lags or simplify the equations

    # pi_t1
    Γ0[eq["pi1"], endo["pi_t1"]] = 1.
    Γ1[eq["pi1"], endo["pi_t"]] = 1.

    # pi_t2
    Γ0[eq["pi2"], endo["pi_t2"]] = 1.
    Γ1[eq["pi2"], endo["pi_t1"]] = 1.

    # pi_a
    Γ0[eq["pi_a"], endo["pi_a_t"]] = 1.
    Γ0[eq["pi_a"], endo["pi_t"]] = -1.
    Γ0[eq["pi_a"], endo["pi_t1"]] = -1.
    Γ0[eq["pi_a"], endo["pi_t2"]] = -1.
    Γ1[eq["pi_a"], endo["pi_t2"]] = 1.

    # Rt1
    Γ0[eq["Rt1"], endo["R_t1"]] = 1.
    Γ1[eq["Rt1"], endo["R_t"]] = 1.

    # E_z
    Γ0[eq["eq_Ez"], endo["E_z"]] = 1.
    Γ0[eq["eq_Ez"], endo["ztil_t"]] = -(Θ.ρ_z-1)/(1-Θ.alp)
    Γ0[eq["eq_Ez"], endo["zp_t"]] = -Θ.ρ_zp



    ### EXOGENOUS SHOCKS ###

    # Neutral technology
    Γ0[eq["eq_z"], endo["z_t"]] = 1.
    Γ1[eq["eq_z"], endo["ztil_t"]] = (Θ.ρ_z - 1)/(1 - Θ.alp)
    Γ0[eq["eq_z"], endo["zp_t"]] = -1.
    Ψ[eq["eq_z"], exo["z_sh"]] = 1/(1 - Θ.alp)

    Γ0[eq["eq_ztil"], endo["ztil_t"]] = 1.
    Γ1[eq["eq_ztil"], endo["ztil_t"]] = Θ.ρ_z
    Ψ[eq["eq_ztil"], exo["z_sh"]] = 1.

    # Long-run changes to productivity
    Γ0[eq["eq_zp"], endo["zp_t"]] = 1.
    Γ1[eq["eq_zp"], endo["zp_t"]] = Θ.ρ_zp
    Ψ[eq["eq_zp"], exo["zp_sh"]] = 1.

    # Government spending
    Γ0[eq["eq_g"], endo["g_t"]] = 1.
    Γ1[eq["eq_g"], endo["g_t"]] = Θ.ρ_g
    Ψ[eq["eq_g"], exo["g_sh"]] = 1.
    Ψ[eq["eq_g"], exo["z_sh"]] = Θ.eta_gz

    # Asset shock
    Γ0[eq["eq_b"], endo["b_t"]] = 1.
    Γ1[eq["eq_b"], endo["b_t"]] = Θ.ρ_b
    Ψ[eq["eq_b"], exo["b_sh"]] = 1.

    # Investment-specific technology
    Γ0[eq["eq_mu"], endo["mu_t"]] = 1.
    Γ1[eq["eq_mu"], endo["mu_t"]] = Θ.ρ_mu
    Ψ[eq["eq_mu"], exo["mu_sh"]] = 1.

    # Price mark-up shock
    Γ0[eq["eq_laf"], endo["laf_t"]] = 1.
    Γ1[eq["eq_laf"], endo["laf_t"]] = Θ.ρ_laf
    Γ1[eq["eq_laf"], endo["laf_t1"]] = -Θ.eta_laf
    Ψ[eq["eq_laf"], exo["laf_sh"]] = 1.

    Γ0[eq["eq_laf1"], endo["laf_t1"]] = 1.
    Ψ[eq["eq_laf1"], exo["laf_sh"]] = 1.

    # Wage mark-up shock
    Γ0[eq["eq_law"], endo["law_t"]] = 1.
    Γ1[eq["eq_law"], endo["law_t"]] = Θ.ρ_law
    Γ1[eq["eq_law"], endo["law_t1"]] = -Θ.eta_law
    Ψ[eq["eq_law"], exo["law_sh"]] = 1.

    Γ0[eq["eq_law1"], endo["law_t1"]] = 1.
    Ψ[eq["eq_law1"], exo["law_sh"]] = 1.

    # Monetary policy shock
    Γ0[eq["eq_rm"], endo["rm_t"]] = 1.
    Γ1[eq["eq_rm"], endo["rm_t"]] = Θ.ρ_rm
    Ψ[eq["eq_rm"], exo["rm_sh"]] = 1.



    ### Financial frictions

    # sigw shock
    Γ0[eq["eq_sigw"], endo["sigw_t"]] = 1.
    Γ1[eq["eq_sigw"], endo["sigw_t"]] = Θ.ρ_sigw
    Ψ[eq["eq_sigw"], exo["sigw_sh"]] = 1.

    # mue shock
    Γ0[eq["eq_mue"], endo["mue_t"]] = 1.
    Γ1[eq["eq_mue"], endo["mue_t"]] = Θ.ρ_mue
    Ψ[eq["eq_mue"], exo["mue_sh"]] = 1.

    # gamm shock
    Γ0[eq["eq_gamm"], endo["gamm_t"]] = 1.
    Γ1[eq["eq_gamm"], endo["gamm_t"]] = Θ.ρ_gamm
    Ψ[eq["eq_gamm"], exo["gamm_sh"]] = 1.

    # Long-term inflation expectations
    Γ0[eq["eq_pist"], endo["pist_t"]] = 1.
    Γ1[eq["eq_pist"], endo["pist_t"]] = Θ.ρ_pist
    Ψ[eq["eq_pist"], exo["pist_sh"]] = 1.

    # Anticipated policy shocks
    if spec["nant"] > 0

        # This section adds the anticipated shocks. There is one state for all the
        # anticipated shocks that will hit in a given period (i.e. rm_tl2 holds those that
        # will hit in two periods), and the equations are set up so that rm_tl2 last period
        # will feed into rm_tl1 this period (and so on for other numbers), and last period's
        # rm_tl1 will feed into the rm_t process (and affect the Taylor Rule this period).

        Γ1[eq["eq_rm"], endo["rm_tl1"]] = 1.
        Γ0[eq["eq_rml1"], endo["rm_tl1"]] = 1.
        Ψ[eq["eq_rml1"], exo["rm_shl1"]] = 1.

        if spec["nant"] > 1
            for i = 2:spec["nant"]
                Γ1[eq["eq_rml$(i-1)"], endo["rm_tl$i"]] = 1.
                Γ0[eq["eq_rml$i"], endo["rm_tl$i"]] = 1.
                Ψ[eq["eq_rml$i"], exo["rm_shl$i"]] = 1.
            end
        end
    end



    ### EXPECTATION ERRORS ###

    ### E(c)

    # Sticky prices and wages
    Γ0[eq["eq_Ec"], endo["c_t"]] = 1.
    Γ1[eq["eq_Ec"], endo["E_c"]] = 1.
    Π[eq["eq_Ec"], ex["Ec_sh"]] = 1.

    # Flexible prices and wages
    Γ0[eq["eq_Ec_f"], endo["c_f_t"]] = 1.
    Γ1[eq["eq_Ec_f"], endo["E_c_f"]] = 1.
    Π[eq["eq_Ec_f"], ex["Ec_f_sh"]] = 1.



    ### E(q)

    # Sticky prices and wages
    Γ0[eq["eq_Eqk"], endo["qk_t"]] = 1.
    Γ1[eq["eq_Eqk"], endo["E_qk"]] = 1.
    Π[eq["eq_Eqk"], ex["Eqk_sh"]] = 1.

    # Flexible prices and wages
    Γ0[eq["eq_Eqk_f"], endo["qk_f_t"]] = 1.
    Γ1[eq["eq_Eqk_f"], endo["E_qk_f"]] = 1.
    Π[eq["eq_Eqk_f"], ex["Eqk_f_sh"]] = 1.



    ### E(i)

    # Sticky prices and wages
    Γ0[eq["eq_Ei"], endo["i_t"]] = 1.
    Γ1[eq["eq_Ei"], endo["E_i"]] = 1.
    Π[eq["eq_Ei"], ex["Ei_sh"]] = 1.

    # Flexible prices and wages
    Γ0[eq["eq_Ei_f"], endo["i_f_t"]] = 1.
    Γ1[eq["eq_Ei_f"], endo["E_i_f"]] = 1.
    Π[eq["eq_Ei_f"], ex["Ei_f_sh"]] = 1.



    ### E(pi)

    # Sticky prices and wages
    Γ0[eq["eq_Epi"], endo["pi_t"]] = 1.
    Γ1[eq["eq_Epi"], endo["E_pi"]] = 1.
    Π[eq["eq_Epi"], ex["Epi_sh"]] = 1.



    ### E(l)

    # Sticky prices and wages
    Γ0[eq["eq_EL"], endo["L_t"]] = 1.
    Γ1[eq["eq_EL"], endo["E_L"]] = 1.
    Π[eq["eq_EL"], ex["EL_sh"]] = 1.

    # Flexible prices and wages
    Γ0[eq["eq_EL_f"], endo["L_f_t"]] = 1.
    Γ1[eq["eq_EL_f"], endo["E_L_f"]] = 1.
    Π[eq["eq_EL_f"], ex["EL_f_sh"]] = 1.



    ### E(rk)

    # Sticky prices and wages
    Γ0[eq["eq_Erk"], endo["rk_t"]] = 1.
    Γ1[eq["eq_Erk"], endo["E_rk"]] = 1.
    Π[eq["eq_Erk"], ex["Erk_sh"]] = 1.

    # Flexible prices and wages
    Γ0[eq["eq_Erk_f"], endo["rk_f_t"]] = 1.
    Γ1[eq["eq_Erk_f"], endo["E_rk_f"]] = 1.
    Π[eq["eq_Erk_f"], ex["Erk_f_sh"]] = 1.



    ### E(w)

    # Sticky prices and wages
    Γ0[eq["eq_Ew"], endo["w_t"]] = 1.
    Γ1[eq["eq_Ew"], endo["E_w"]] = 1.
    Π[eq["eq_Ew"], ex["Ew_sh"]] = 1.



    ### E(Rktil)

    # Sticky prices and wages
    Γ0[eq["eq_ERktil"], endo["Rktil_t"]] = 1.
    Γ1[eq["eq_ERktil"], endo["E_Rktil"]] = 1.
    Π[eq["eq_ERktil"], ex["ERktil_sh"]] = 1.



    return Γ0, Γ1, C, Ψ, Π
end





# Assign measurement equation : X_t = ZZ*S_t + DD + u_t
# where u_t = eta_t+MM* eps_t with var(eta_t) = EE
# where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
function measurement(model::Model990, TTT::Matrix, RRR::Matrix, CCC::Matrix; shocks::Bool = true)
    spec = model.spec
    Θ = model.Θ
    ind = model.ind

    endo = ind.endostates
    exo = ind.exoshocks
    obs = ind.observables
    
    # If shocks = true, then return measurement equation matrices with rows and columns for anticipated policy shocks
    if shocks
        n_observables = spec["n_observables"]
        n_states = spec["n_states_aug"]
        n_exoshocks = spec["n_exoshocks"]
        endo_addl = ind.endostates_postgensys
    else
        nant = spec["nant"]
        n_observables = spec["n_observables"] - nant
        n_states = spec["n_states_aug"] - nant
        n_exoshocks = spec["n_exoshocks"] - nant
        endo_addl = [key => ind.endostates_postgensys[key] - nant for key in keys(ind.endostates_postgensys)]
    end

    ZZ = zeros(n_observables, n_states)
    DD = zeros(n_observables, 1)
    MM = zeros(n_observables, n_exoshocks)
    EE = zeros(n_observables, n_observables)
    QQ =  zeros(n_exoshocks, n_exoshocks)



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
    if shocks
        for i = 1:spec["nant"]
            ZZ[obs["R_n$i"], :] = ZZ[obs["R_n"], :]*(TTT^i)
            DD[obs["R_n$i"]] = Θ.Rstarn
            QQ[exo["rm_shl$i"], exo["rm_shl$i"]] = Θ.(parse("σ_rm$i"))^2
        end
    end

    return ZZ, DD, QQ, EE, MM
end
