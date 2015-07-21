# Return a ModelInds object, which has four fields: each is a dictionary mapping states or
#   equations to model indices
# Updates spec with n_states, etc.
function ModelInds(spec::Dict{String, Any})
    # Endogenous states
    endostates = ["y_t", "c_t", "i_t", "qk_t", "k_t", "kbar_t", "u_t", "rk_t", "Rktil_t", "n_t", "mc_t", "pi_t", "muw_t", "w_t", "L_t", "R_t", "g_t", "b_t", "mu_t", "z_t", "laf_t", "laf_t1", "law_t", "law_t1", "rm_t", "sigw_t", "mue_t", "gamm_t", "pist_t", "E_c", "E_qk", "E_i", "E_pi", "E_L", "E_rk", "E_w", "E_Rktil", "y_f_t", "c_f_t", "i_f_t", "qk_f_t", "k_f_t", "kbar_f_t", "u_f_t", "rk_f_t", "w_f_t", "L_f_t", "r_f_t", "E_c_f", "E_qk_f", "E_i_f", "E_L_f", "E_rk_f", "ztil_t", "pi_t1", "pi_t2", "pi_a_t", "R_t1", "zp_t", "E_z", ["rm_tl$i" for i = 1:spec["nant"]]]

    # Exogenous shocks
    exoshocks = ["g_sh", "b_sh", "mu_sh", "z_sh", "laf_sh", "law_sh", "rm_sh", "sigw_sh", "mue_sh", "gamm_sh", "pist_sh", "lr_sh", "zp_sh", "tfp_sh", "gdpdef_sh", "pce_sh", ["rm_shl$i" for i = 1:spec["nant"]]]

    # Expectations shocks
    expshocks = ["Ec_sh", "Eqk_sh", "Ei_sh", "Epi_sh", "EL_sh", "Erk_sh", "Ew_sh", "ERktil_sh", "Ec_f_sh", "Eqk_f_sh", "Ei_f_sh", "EL_f_sh", "Erk_f_sh"]

    # Equilibrium conditions
    eqconds = ["euler", "inv", "capval", "spread", "nevol", "output", "caputl", "capsrv", "capev", "mkupp", "phlps", "caprnt", "msub", "wage", "mp", "res", "eq_g", "eq_b", "eq_mu", "eq_z", "eq_laf", "eq_law", "eq_rm", "eq_sigw", "eq_mue", "eq_gamm", "eq_laf1", "eq_law1", "eq_Ec", "eq_Eqk", "eq_Ei", "eq_Epi", "eq_EL", "eq_Erk", "eq_Ew", "eq_ERktil", "euler_f", "inv_f", "capval_f", "output_f", "caputl_f", "capsrv_f", "capev_f", "mkupp_f", "caprnt_f", "msub_f", "res_f", "eq_Ec_f", "eq_Eqk_f", "eq_Ei_f", "eq_EL_f", "eq_Erk_f", "eq_ztil", "eq_pist", "pi1", "pi2", "pi_a", "Rt1", "eq_zp", "eq_Ez", ["eq_rml$i" for i=1:spec["nant"]]]
    
    # Additional states added after solving model
    # Lagged states and observables measurement error
    endostates_postgensys = ["y_t1", "c_t1", "i_t1", "w_t1", "pi_t1", "L_t1", "Et_pi_t", "lr_t", "tfp_t", "e_gdpdef", "e_pce", "u_t1"]

    # Measurement equation observables
    observables = ["g_y", # quarterly output growth
                   "hoursg", # aggregate hours growth
                   "g_w", # real wage growth
                   "pi_gdpdef", # inflation (GDP deflator)
                   "pi_pce", # inflation (core PCE)
                   "R_n", # nominal interest rate
                   "g_c", # consumption growth
                   "g_i", # investment growth
                   "sprd", # spreads
                   "pi_long", # 10-year inflation expectation
                   "R_long", # long-term rate
                   "tfp", # total factor productivity
                   ["R_n$i" for i = 1:spec["nant"]]] # compounded nominal rates

    # Update spec dictionary with number of states, etc.
    spec["n_states"] = length(endostates)
    spec["n_states_aug"] = length(endostates) + length(endostates_postgensys)
    spec["n_exoshocks"] = length(exoshocks)
    spec["n_expshocks"] = length(expshocks)
    spec["n_eqconds"] = length(eqconds)
    spec["n_observables"] = length(observables)

    return ModelInds(makedict(endostates), makedict(exoshocks), makedict(expshocks),
                     makedict(eqconds), makedict(endostates_postgensys; start=length(endostates)),
                     makedict(observables))
end
