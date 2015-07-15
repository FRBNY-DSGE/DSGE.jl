using ..AbstractModel
import ..AbstractModel: ModelInds

# Return a ModelInds object, which has four fields: each is a dictionary mapping states or
#   equations to model indices
function AbstractModel.ModelInds(spec::Dict{String, Any})
    endostates = ["y_t",
                  "c_t",
                  "i_t",
                  "qk_t",
                  "k_t",
                  "kbar_t",
                  "u_t",
                  "rk_t",
                  "Rktil_t",
                  "n_t",
                  "mc_t",
                  "pi_t",
                  "muw_t", # exo
                  "w_t",
                  "L_t",
                  "R_t",
                  "g_t", # exo
                  "b_t", # exo
                  "mu_t", # exo
                  "z_t", # exo
                  "laf_t", # exo
                  "laf_t1", # exo
                  "law_t", # exo
                  "law_t1", # exo
                  "rm_t", # exo
                  "sigw_t", # exo - FF
                  "mue_t", # exo - FF
                  "gamm_t", # exo - FF
                  "pist_t", # exo
                  "E_c", # exp
                  "E_qk", # exp
                  "E_i", # exp
                  "E_pi", # exp
                  "E_L", # exp
                  "E_rk", # exp
                  "E_w", # exp
                  "E_Rktil", # exp - FF
                  "y_f_t",
                  "c_f_t",
                  "i_f_t",
                  "qk_f_t",
                  "k_f_t",
                  "kbar_f_t",
                  "u_f_t",
                  "rk_f_t",
                  "w_f_t",
                  "L_f_t",
                  "r_f_t",
                  "E_c_f", # exp
                  "E_qk_f", # exp
                  "E_i_f", # exp
                  "E_L_f", # exp
                  "E_rk_f", # exp
                  "ztil_t",
                  "pi_t1",
                  "pi_t2",
                  "pi_a_t",
                  "R_t1",
                  "zp_t",
                  "E_z",
                  ["rm_tl$i" for i = 1:spec["nant"]]]
    
    exoshocks = ["g_sh",
                 "b_sh",
                 "mu_sh",
                 "z_sh",
                 "laf_sh",
                 "law_sh",
                 "rm_sh",
                 "sigw_sh", # FF
                 "mue_sh", # FF
                 "gamm_sh", # FF
                 "pist_sh",
                 "lr_sh",
                 "zp_sh",
                 "tfp_sh",
                 "gdpdef_sh",
                 "pce_sh",
                 ["rm_shl$i" for i = 1:spec["nant"]]]

    expshocks = ["Ec_sh",
                 "Eqk_sh",
                 "Ei_sh",
                 "Epi_sh",
                 "EL_sh",
                 "Erk_sh",
                 "Ew_sh",
                 "ERktil_sh",
                 "Ec_f_sh",
                 "Eqk_f_sh",
                 "Ei_f_sh",
                 "EL_f_sh",
                 "Erk_f_sh"]

    equations = ["euler",
                 "inv",
                 "capval",
                 "spread",
                 "nevol",
                 "output",
                 "caputl",
                 "capsrv",
                 "capev",
                 "mkupp",
                 "phlps",
                 "caprnt",
                 "msub",
                 "wage",
                 "mp",
                 "res",
                 "eq_g",
                 "eq_b",
                 "eq_mu",
                 "eq_z",
                 "eq_laf",
                 "eq_law",
                 "eq_rm",
                 "eq_sigw",
                 "eq_mue",
                 "eq_gamm",
                 "eq_laf1",
                 "eq_law1",
                 "eq_Ec",
                 "eq_Eqk",
                 "eq_Ei",
                 "eq_Epi",
                 "eq_EL",
                 "eq_Erk",
                 "eq_Ew",
                 "eq_ERktil",
                 "euler_f",
                 "inv_f",
                 "capval_f",
                 "output_f",
                 "caputl_f",
                 "capsrv_f",
                 "capev_f",
                 "mkupp_f",
                 "caprnt_f",
                 "msub_f",
                 "res_f",
                 "eq_Ec_f",
                 "eq_Eqk_f",
                 "eq_Ei_f",
                 "eq_EL_f",
                 "eq_Erk_f",
                 "eq_ztil",
                 "eq_pist",
                 "pi1",
                 "pi2",
                 "pi_a",
                 "Rt1",
                 "eq_zp",
                 "eq_Ez",
                 ["eq_rml$i" for i=1:spec["nant"]]]

    # TODO: incorporate measurement equation observables
    observables = ["g_y", # output growth
                   "g_l", # aggregate hours growth
                   "g_w", # real wage growth
                   "gdpdef", # GDP deflator
                   "pi_pce", # core PCE inflation
                   "R", # interest rate
                   "g_c", # consumption growth
                   "g_i", # investment growth
                   "sprd", # spread
                   "pi_long", # long-term inflation
                   "R_long", # long-term rate
                   "tfp"] # total factor productivity

    return ModelInds(makedict(endostates), makedict(exoshocks), makedict(expshocks),
                     makedict(equations))
end
