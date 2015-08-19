# The given fields define the entire model structure.
# We can then concisely pass around a Model object to the remaining steps of the model
#   (solve, estimate, and forecast).
type Model990 <: AbstractDSGEModel
    parameters::Vector{Param}                       # vector of all of the model parameters
    parameters_fixed::Vector{Param}                 # vector of all "permanently fixed" model parameters
    steady_state::Vector                            # model steady-state values
    keys::Dict{Symbol,Int}                          # human-readable names for all the model
                                                    # parameters and steady-num_states

    endogenous_states::Dict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::Dict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::Dict{Symbol,Int}               #
    equilibrium_conditions::Dict{Symbol,Int}        #
    endogenous_states_postgensys::Dict{Symbol,Int}  #
    observables::Dict{Symbol,Int}                   #

    num_anticipated_shocks::Int                     # Number of anticipated policy shocks
    num_anticipated_shocks_padding::Int             # Padding for nant
    num_anticipated_lags::Int                       # Number of periods back to incorporate zero bound expectations
    num_presample_periods::Int                      # Number of periods in the presample

    reoptimize::Bool                                # Reoptimize the posterior mode
    recalculate_hessian::Bool                       # Recalculate the hessian at the mode
    num_mh_simulations::Int                         # 
    num_mh_blocks::Int                              #
    num_mh_burn::Int                                #
    mh_thinning_step::Int                           #
end

description(m::Model990) = "This is some model that we're trying to make work."

function initialise_model_indices!(m::Model990)
    # Endogenous states
    endogenous_states = [
      :y_t, :c_t, :i_t, :qk_t, :k_t, :kbar_t, :u_t, :rk_t, :Rktil_t, :n_t, :mc_t,
      :pi_t, :muw_t, :w_t, :L_t, :R_t, :g_t, :b_t, :mu_t, :z_t, :laf_t, :laf_t1,
      :law_t, :law_t1, :rm_t, :sigw_t, :mue_t, :gamm_t, :pist_t, :E_c, :E_qk, :E_i,
      :E_pi, :E_L, :E_rk, :E_w, :E_Rktil, :y_f_t, :c_f_t, :i_f_t, :qk_f_t, :k_f_t,
      :kbar_f_t, :u_f_t, :rk_f_t, :w_f_t, :L_f_t, :r_f_t, :E_c_f, :E_qk_f, :E_i_f,
      :E_L_f, :E_rk_f, :ztil_t, :pi_t1, :pi_t2, :pi_a_t, :R_t1, :zp_t, :E_z,
      [symbol("rm_tl$i") for i = 1:num_anticipated_shocks(m)]]

    # Exogenous shocks
    exogenous_shocks = [
      :g_sh, :b_sh, :mu_sh, :z_sh, :laf_sh, :law_sh, :rm_sh, :sigw_sh, :mue_sh,
      :gamm_sh, :pist_sh, :lr_sh, :zp_sh, :tfp_sh, :gdpdef_sh, :pce_sh,
      [symbol("rm_shl$i") for i = 1:num_anticipated_shocks(m)]]

    # Expectations shocks
    expected_shocks = [
      : Ec_sh, :Eqk_sh, :Ei_sh, :Epi_sh, :EL_sh, :Erk_sh, :Ew_sh, :ERktil_sh, :Ec_f_sh,
      :Eqk_f_sh, :Ei_f_sh, :EL_f_sh, :Erk_f_sh]

    # Equilibrium conditions
    equilibrium_conditions = [
      :euler, :inv, :capval, :spread, :nevol, :output, :caputl, :capsrv, :capev,
      :mkupp, :phlps, :caprnt, :msub, :wage, :mp, :res, :eq_g, :eq_b, :eq_mu, :eq_z,
      :eq_laf, :eq_law, :eq_rm, :eq_sigw, :eq_mue, :eq_gamm, :eq_laf1, :eq_law1, :eq_Ec,
      :eq_Eqk, :eq_Ei, :eq_Epi, :eq_EL, :eq_Erk, :eq_Ew, :eq_ERktil, :euler_f, :inv_f,
      :capval_f, :output_f, :caputl_f, :capsrv_f, :capev_f, :mkupp_f, :caprnt_f, :msub_f,
      :res_f, :eq_Ec_f, :eq_Eqk_f, :eq_Ei_f, :eq_EL_f, :eq_Erk_f, :eq_ztil, :eq_pist,
      :pi1, :pi2, :pi_a, :Rt1, :eq_zp, :eq_Ez, [symbol("eq_rml$i") for i=1:num_anticipated_shocks(m)]]

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_postgensys = [
      :y_t1, :c_t1, :i_t1, :w_t1, :pi_t1, :L_t1, :Et_pi_t, :lr_t, :tfp_t, :e_gdpdef,
      :e_pce, :u_t1]

    # Measurement equation observables
    observables = [:g_y,         # quarterly output growth
                   :hoursg,      # aggregate hours growth
                   :g_w,         # real wage growth
                   :pi_gdpdef,   # inflation (GDP deflator)
                   :pi_pce,      # inflation (core PCE)
                   :R_n,         # nominal interest rate
                   :g_c,         # consumption growth
                   :g_i,         # investment growth
                   :sprd,        # spreads
                   :pi_long,     # 10-year inflation expectation
                   :R_long,      # long-term rate
                   :tfp,         # total factor productivity
                   [symbol("R_n$i") for i=1:num_anticipated_shocks(m)]] # compounded nominal rates

    for (i,k) in enumerate(endogenous_states);            m.endogenous_states[k]            = i end
    for (i,k) in enumerate(exogenous_shocks);             m.exogenous_shocks[k]             = i end
    for (i,k) in enumerate(expected_shocks);              m.expected_shocks[k]              = i end
    for (i,k) in enumerate(equilibrium_conditions);       m.equilibrium_conditions[k]       = i end
    for (i,k) in enumerate(endogenous_states);            m.endogenous_states[k]            = i end
    for (i,k) in enumerate(endogenous_states_postgensys); m.endogenous_states_postgensys[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                  m.observables[k]                  = i end
end

function initialise_model_parameters!(m::Model990)
    m[:alp     ] = Param(0.1596, false, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999))
    m[:zeta_p  ] = Param(0.8940, false, (1e-5, 0.999), BetaAlt(0.5, 0.1), 1, (1e-5, 0.999))
    m[:iota_p  ] = Param(0.1865, false, (1e-5, 0.999), BetaAlt(0.5, 0.15), 1, (1e-5, 0.999))
    m[:del     ] = Param(0.025) # omit from parameter vector
    m[:ups     ] = Param(1.000, true, (0., 10.), GammaAlt(1., 0.5), 2, (1e-5, 0.))
    m[:Bigphi  ] = Param(1.1066, false, (1., 10.), Normal(1.25, 0.12), 2, (1.00, 10.00))
    m[:s2      ] = Param(2.7314, false, (-15., 15.), Normal(4., 1.5), 0, (-15., 15.))
    m[:h       ] = Param(0.5347, false, (1e-5, 0.999), BetaAlt(0.7, 0.1), 1, (1e-5, 0.999))
    m[:ppsi    ] = Param(0.6862, false, (1e-5, 0.999), BetaAlt(0.5, 0.15), 1, (1e-5, 0.999))
    m[:nu_l    ] = Param(2.5975, false, (1e-5, 10.), Normal(2, 0.75), 2, (1e-5, 10.))
    m[:zeta_w  ] = Param(0.9291, false, (1e-5, 0.999), BetaAlt(0.5, 0.1), 1, (1e-5, 0.999))
    m[:iota_w  ] = Param(0.2992, false, (1e-5, 0.999), BetaAlt(0.5, 0.15), 1, (1e-5, 0.999))
    m[:law     ] = Param(1.5) # omit from parameter vector

    m[:bet     ] = Param(0.1402, scalefunction = x -> 1/(1 + x/100), false, (1e-5, 10.), GammaAlt(0.25, 0.1), 2, (1e-5, 10.))
    m[:psi1    ] = Param(1.3679, false, (1e-5, 10.), Normal(1.5, 0.25), 2, (1e-5, 10.00))
    m[:psi2    ] = Param(0.0388, false, (-0.5, 0.5), Normal(0.12, 0.05), 0, (-0.5, 0.5))
    m[:psi3    ] = Param(0.2464, false, (-0.5, 0.5), Normal(0.12, 0.05), 0, (-0.5, 0.5))
    m[:pistar  ] = Param(0.5000, scalefunction = x -> 1 + x/100, true, (1e-5, 10.), GammaAlt(0.75, 0.4), 2, (1e-5, 10.))
    m[:sigmac  ] = Param(0.8719, false, (1e-5, 10.), Normal(1.5, 0.37), 2, (1e-5, 10.))
    m[:rho     ] = Param(0.7126, false, (1e-5, 0.999), BetaAlt(0.75, 0.10), 1, (1e-5, 0.999))
    m[:epsp    ] = Param(10.) # omit from parameter vector
    m[:epsw    ] = Param(10.) # omit from parameter vector

    # financial frictions parameters
    m[:Fom     ] = Param(0.0300, scalefunction = x -> 1 - (1-x)^0.25, true, (1e-5, 0.99999), BetaAlt(0.03, 0.01), 1, (1e-5, 0.99))
    m[:sprd    ] = Param(1.7444, scalefunction = x -> (1 + x/100)^0.25, false, (0., 100.), GammaAlt(2., 0.1), 2, (1e-5, 0.))
    m[:zeta_spb] = Param(0.0559, false, (1e-5, 0.99999), BetaAlt(0.05, 0.005), 1, (1e-5, 0.99))
    m[:gammstar] = Param(0.9900, true, (1e-5, 0.99999), BetaAlt(0.99, 0.002), 1, (1e-5, 0.99))

    # exogenous processes - level
    m[:gam     ] = Param(0.3673, scalefunction = x -> x/100, false, (-5., 5.), Normal(0.4, 0.1), 0, (-5.0, 5.0))
    m[:Lmean   ] = Param(-45.9364, false, (-1000., 1000.), Normal(-45, 5), 0, (-1000., 1000.))
    m[:gstar   ] = Param(0.18) # omit from parameter vector

    # exogenous processes - autocorrelation
    m[:ρ_g     ] = Param(0.9863, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_b     ] = Param(0.9410, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_mu    ] = Param(0.8735, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_z     ] = Param(0.9446, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_laf   ] = Param(0.8827, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_law   ] = Param(0.3884, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_rm    ] = Param(0.2135, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_sigw  ] = Param(0.9898, false, (1e-5, 0.99999), BetaAlt(0.75, 0.15), 1, (1e-5, 0.99))
    m[:ρ_mue   ] = Param(0.7500, true, (1e-5, 0.99999), BetaAlt(0.75, 0.15), 1, (1e-5, 0.99))
    m[:ρ_gamm  ] = Param(0.7500, true, (1e-5, 0.99999), BetaAlt(0.75, 0.15), 1, (1e-5, 0.99))
    m[:ρ_pist  ] = Param(0.9900, true, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_lr    ] = Param(0.6936, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_zp    ] = Param(0.8910, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_tfp   ] = Param(0.1953, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_gdpdef] = Param(0.5379, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))
    m[:ρ_pce   ] = Param(0.2320, false, (1e-5, 0.999), BetaAlt(0.5, 0.2), 1, (1e-5, 0.999))

    # exogenous processes - standard deviation
    m[:σ_g     ] = Param(2.5230, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_b     ] = Param(0.0292, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_mu    ] = Param(0.4559, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_z     ] = Param(0.6742, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_laf   ] = Param(0.1314, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_law   ] = Param(0.3864, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_rm    ] = Param(0.2380, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_sigw  ] = Param(0.0428, false, (1e-7, 100.), RootInverseGamma(4.00, 0.05), 2, (1e-5, 0.))
    m[:σ_mue   ] = Param(0., true, (1e-7, 100.), RootInverseGamma(4.00, 0.05), 2, (1e-5, 0.))
    m[:σ_gamm  ] = Param(0., true, (1e-7, 100.), RootInverseGamma(4.00, 0.01), 2, (1e-5, 0.))
    m[:σ_pist  ] = Param(0.0269, false, (1e-8, 5.), RootInverseGamma(6., 0.03), 2, (1e-8, 5.))
    m[:σ_lr    ] = Param(0.1766, false, (1e-8, 10.), RootInverseGamma(2., 0.75), 2, (1e-8, 5.))
    m[:σ_zp    ] = Param(0.1662, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_tfp   ] = Param(0.9391, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_gdpdef] = Param(0.1575, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    m[:σ_pce   ] = Param(0.0999, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))

    # standard deviations of the anticipated policy shocks
    for i = 1:num_anticipated_shocks_padding(m)
        m[symbol("σ_rm$i")] = if i < 13
            Param(0.20, false, (1e-7, 100.), RootInverseGamma(4.00, 0.2), 2, (1e-5, 0.))
        else
            Param(0.0, true, (1e-7, 100.), RootInverseGamma(4.00, 0.2), 2, (1e-5, 0.))
        end
    end

    m[:eta_gz      ] = Param(0.8400, false, (1e-5, 0.999), BetaAlt(0.50, 0.20), 1, (1e-5, 0.999))
    m[:eta_laf     ] = Param(0.7892, false, (1e-5, 0.999), BetaAlt(0.50, 0.20), 1, (1e-5, 0.999))
    m[:eta_law     ] = Param(0.4226, false, (1e-5, 0.999), BetaAlt(0.50, 0.20), 1, (1e-5, 0.999))

    m[:modelalp_ind] = Param(0.0000, true, (0.000, 1.000), BetaAlt(0.50, 0.20), 0, (0., 0.))
    m[:gamm_gdpdef ] = Param(1.0354, false, (-10., 10.), Normal(1.00, 2.), 0, (-10., -10.))
    m[:del_gdpdef  ] = Param(0.0181, false, (-9.1, 9.1), Normal(0.00, 2.), 0, (-10., -10.))
end

function Model990()
    parameter_keys, parameter_fixed_keys, steady_state_keys = (
            # parameter keys
            [:alp, :zeta_p, :iota_p, :ups, :Bigphi, :s2, :h, :ppsi, :nu_l, :zeta_w, :iota_w,
            :bet, :psi1, :psi2, :psi3, :pistar, :sigmac, :rho, :Fom, :sprd,
            :zeta_spb, :gammstar, :gam, :Lmean, :ρ_g, :ρ_b, :ρ_mu, :ρ_z, :ρ_laf, :ρ_law,
            :ρ_rm, :ρ_sigw, :ρ_mue, :ρ_gamm, :ρ_pist, :ρ_lr, :ρ_zp, :ρ_tfp, :ρ_gdpdef, :ρ_pce,
            :σ_g, :σ_b, :σ_mu, :σ_z, :σ_laf, :σ_law, :σ_rm, :σ_sigw, :σ_mue, :σ_gamm, :σ_pist,
            :σ_lr, :σ_zp, :σ_tfp, :σ_gdpdef, :σ_pce, :σ_rm1, :σ_rm2, :σ_rm3, :σ_rm4, :σ_rm5,
            :σ_rm6, :σ_rm7, :σ_rm8, :σ_rm9, :σ_rm10, :σ_rm11, :σ_rm12, :σ_rm13, :σ_rm14, :σ_rm15,
            :σ_rm16, :σ_rm17, :σ_rm18, :σ_rm19, :σ_rm20, :eta_gz, :eta_laf, :eta_law, :modelalp_ind,
            :gamm_gdpdef, :del_gdpdef],

            # parameter fixed keys
            [:del, :law, :epsp, :epsw, :gstar],

            # steady state keys
            [:zstar, :rstar, :Rstarn, :rkstar, :wstar, :Lstar, :kstar, :kbarstar, :istar, :ystar,
            :cstar, :wl_c, :nstar, :vstar, :zeta_spsigw, :zeta_spmue, :zeta_nRk, :zeta_nR, :zeta_nqk,
            :zeta_nn, :zeta_nmue, :zeta_nsigw]
            )

    # initialise human-readable keys for variablesiables
    keylist = Dict{Symbol,Int}()
    for (i,k) in enumerate(vcat(parameter_keys,parameter_fixed_keys,steady_state_keys))
        keylist[k] = i
    end

    # initialise vector to store actual values
    parameters       = @compat Vector{Param}(length(parameter_keys))
    parameters_fixed = @compat Vector{Param}(length(parameter_fixed_keys))
    steady_state     = @compat Vector{Any}(length(steady_state_keys))

    # Model-specific specifications
    num_anticipated_shocks          = 6
    num_anticipated_shocks_padding  = 20
    num_anticipated_lags            = 24
    num_presample_periods           = 2 # TODO: This should be set when the data are read in

    # Estimation specifications
    reoptimize                      = false
    recalculate_hessian             = false
    num_mh_simulations              = 10000
    num_mh_blocks                   = 22
    num_mh_burn                     = 2
    mh_thinning_step                = 5
    num_mh_simulations_test         = 10
    num_mh_blocks_test              = 1
    num_mh_burn_test                = 0

    # initialise empty model
    m = Model990(
            parameters,
            parameters_fixed,
            steady_state,
            keylist,
            Dict{Symbol,Int}(),
            Dict{Symbol,Int}(),
            Dict{Symbol,Int}(),
            Dict{Symbol,Int}(),
            Dict{Symbol,Int}(),
            Dict{Symbol,Int}(),

            num_anticipated_shocks,
            num_anticipated_shocks_padding,
            num_anticipated_lags,
            num_presample_periods,

            reoptimize,
            recalculate_hessian,
            num_mh_simulations,
            num_mh_blocks,
            num_mh_burn,
            mh_thinning_step)

    initialise_model_parameters!(m)
    initialise_model_indices!(m)

    return steadystate!(m)
end

# functions that are used to compute financial frictions
# steady-state values from parameter values
@inline function ζ_spb_fn(z, σ, sprd)
    zetaratio = ζ_bω_fn(z, σ, sprd)/ζ_zω_fn(z, σ, sprd)
    nk = nk_fn(z, σ, sprd)
    return -zetaratio/(1-zetaratio)*nk/(1-nk)
end

@inline function ζ_bω_fn(z, σ, sprd)
    nk          = nk_fn(z, σ, sprd)
    μstar       = μ_fn(z, σ, sprd)
    ωstar       = ω_fn(z, σ)
    Γstar       = Γ_fn(z, σ)
    Gstar       = G_fn(z, σ)
    dΓ_dωstar   = dΓ_dω_fn(z)
    dG_dωstar   = dG_dω_fn(z, σ)
    d2Γ_dω2star = d2Γ_dω2_fn(z, σ)
    d2G_dω2star = d2G_dω2_fn(z, σ)
    return ωstar*μstar*nk*(d2Γ_dω2star*dG_dωstar - d2G_dω2star*dΓ_dωstar)/
        (dΓ_dωstar - μstar*dG_dωstar)^2/sprd/(1 - Γstar + dΓ_dωstar*(Γstar - μstar*Gstar)/
            (dΓ_dωstar - μstar*dG_dωstar))
end

@inline function ζ_zω_fn(z, σ, sprd)
    μstar = μ_fn(z, σ, sprd)
    return ω_fn(z, σ)*(dΓ_dω_fn(z) - μstar*dG_dω_fn(z, σ))/
        (Γ_fn(z, σ) - μstar*G_fn(z, σ))
end

nk_fn(z, σ, sprd) = 1 - (Γ_fn(z, σ) - μ_fn(z, σ, sprd)*G_fn(z, σ))*sprd
μ_fn(z, σ, sprd)  =
    (1 - 1/sprd)/(dG_dω_fn(z, σ)/dΓ_dω_fn(z)*(1 - Γ_fn(z, σ)) + G_fn(z, σ))
ω_fn(z, σ)        = exp(σ*z - σ^2/2)
G_fn(z, σ)        = cdf(Normal(), z-σ)
Γ_fn(z, σ)        = ω_fn(z, σ)*(1 - cdf(Normal(), z)) + cdf(Normal(), z-σ)
dG_dω_fn(z, σ)    = pdf(Normal(), z)/σ
d2G_dω2_fn(z, σ)  = -z*pdf(Normal(), z)/ω_fn(z, σ)/σ^2
dΓ_dω_fn(z)       = 1 - cdf(Normal(), z)
d2Γ_dω2_fn(z, σ)  = -pdf(Normal(), z)/ω_fn(z, σ)/σ
dG_dσ_fn(z, σ)    = -z*pdf(Normal(), z-σ)/σ
d2G_dωdσ_fn(z, σ) = -pdf(Normal(), z)*(1 - z*(z-σ))/σ^2
dΓ_dσ_fn(z, σ)    = -pdf(Normal(), z-σ)
d2Γ_dωdσ_fn(z, σ) = (z/σ-1)*pdf(Normal(), z)

# (Re)calculates steady-state values
function steadystate!(m::Model990)
    m[:zstar]    = log(1+m[:gam]) + m[:alp]/(1-m[:alp])*log(m[:ups])
    m[:rstar]    = exp(m[:sigmac]*m[:zstar]) / m[:bet]
    m[:Rstarn]   = 100*(m[:rstar]*m[:pistar] - 1)
    m[:rkstar]   = m[:sprd]*m[:rstar]*m[:ups] - (1-m[:del])
    m[:wstar]    = (m[:alp]^m[:alp] * (1-m[:alp])^(1-m[:alp]) * m[:rkstar]^(-m[:alp]) / m[:Bigphi])^(1/(1-m[:alp]))
    m[:Lstar]    = 1.
    m[:kstar]    = (m[:alp]/(1-m[:alp])) * m[:wstar] * m[:Lstar] / m[:rkstar]
    m[:kbarstar] = m[:kstar] * (1+m[:gam]) * m[:ups]^(1 / (1-m[:alp]))
    m[:istar]    = m[:kbarstar] * (1-((1-m[:del])/((1+m[:gam]) * m[:ups]^(1/(1-m[:alp])))))
    m[:ystar]    = (m[:kstar]^m[:alp]) * (m[:Lstar]^(1-m[:alp])) / m[:Bigphi]
    m[:cstar]    = (1-m[:gstar])*m[:ystar] - m[:istar]
    m[:wl_c]     = (m[:wstar]*m[:Lstar])/(m[:cstar]*m[:law])

    # FINANCIAL FRICTIONS ADDITIONS
    # solve for sigmaomegastar and zomegastar
    zwstar = quantile(Normal(), m[:Fom].scaledvalue)
    sigwstar = fzero(sigma -> ζ_spb_fn(zwstar, sigma, m[:sprd]) - m[:zeta_spb], 0.5)

    # evaluate omegabarstar
    omegabarstar = ω_fn(zwstar, sigwstar)

    # evaluate all BGG function elasticities
    Gstar                   = G_fn(zwstar, sigwstar)
    Gammastar               = Γ_fn(zwstar, sigwstar)
    dGdomegastar            = dG_dω_fn(zwstar, sigwstar)
    d2Gdomega2star          = d2G_dω2_fn(zwstar, sigwstar)
    dGammadomegastar        = dΓ_dω_fn(zwstar)
    d2Gammadomega2star      = d2Γ_dω2_fn(zwstar, sigwstar)
    dGdsigmastar            = dG_dσ_fn(zwstar, sigwstar)
    d2Gdomegadsigmastar     = d2G_dωdσ_fn(zwstar, sigwstar)
    dGammadsigmastar        = dΓ_dσ_fn(zwstar, sigwstar)
    d2Gammadomegadsigmastar = d2Γ_dωdσ_fn(zwstar, sigwstar)

    # evaluate mu, nk, and Rhostar
    muestar       = μ_fn(zwstar, sigwstar, m[:sprd])
    nkstar        = nk_fn(zwstar, sigwstar, m[:sprd])
    Rhostar       = 1/nkstar - 1

    # evaluate wekstar and vkstar
    wekstar       = (1-m[:gammstar]/m[:bet])*nkstar - m[:gammstar]/m[:bet]*(m[:sprd]*(1-muestar*Gstar) - 1)
    vkstar        = (nkstar-wekstar)/m[:gammstar]

    # evaluate nstar and vstar
    m[:nstar]       = nkstar*m[:kstar]
    m[:vstar]       = vkstar*m[:kstar]

    # a couple of combinations
    GammamuG      = Gammastar - muestar*Gstar
    GammamuGprime = dGammadomegastar - muestar*dGdomegastar

    # elasticities wrt omegabar
    zeta_bw       = ζ_bω_fn(zwstar, sigwstar, m[:sprd])
    zeta_zw       = ζ_zω_fn(zwstar, sigwstar, m[:sprd])
    zeta_bw_zw    = zeta_bw/zeta_zw

    # elasticities wrt sigw
    zeta_bsigw      = sigwstar * (((1 - muestar*dGdsigmastar/dGammadsigmastar) /
        (1 - muestar*dGdomegastar/dGammadomegastar) - 1)*dGammadsigmastar*m[:sprd] + muestar*nkstar*
            (dGdomegastar*d2Gammadomegadsigmastar - dGammadomegastar*d2Gdomegadsigmastar)/GammamuGprime^2) /
                ((1 - Gammastar)*m[:sprd] + dGammadomegastar/GammamuGprime*(1-nkstar))
    zeta_zsigw      = sigwstar * (dGammadsigmastar - muestar*dGdsigmastar) / GammamuG
    m[:zeta_spsigw] = (zeta_bw_zw*zeta_zsigw - zeta_bsigw) / (1-zeta_bw_zw)

    # elasticities wrt mue
    zeta_bmue      = muestar * (nkstar*dGammadomegastar*dGdomegastar/GammamuGprime+dGammadomegastar*Gstar*m[:sprd]) /
        ((1-Gammastar)*GammamuGprime*m[:sprd] + dGammadomegastar*(1-nkstar))
    zeta_zmue      = -muestar*Gstar/GammamuG
    m[:zeta_spmue] = (zeta_bw_zw*zeta_zmue - zeta_bmue) / (1-zeta_bw_zw)

    # some ratios/elasticities
    Rkstar        = m[:sprd]*m[:pistar]*m[:rstar] # (rkstar+1-delta)/ups*pistar
    zeta_gw       = dGdomegastar/Gstar*omegabarstar
    zeta_Gsigw    = dGdsigmastar/Gstar*sigwstar

    # elasticities for the net worth evolution
    m[:zeta_nRk]    = m[:gammstar]*Rkstar/m[:pistar]/exp(m[:zstar])*(1+Rhostar)*(1 - muestar*Gstar*(1 - zeta_gw/zeta_zw))
    m[:zeta_nR]     = m[:gammstar]/m[:bet]*(1+Rhostar)*(1 - nkstar + muestar*Gstar*m[:sprd]*zeta_gw/zeta_zw)
    m[:zeta_nqk]    = m[:gammstar]*Rkstar/m[:pistar]/exp(m[:zstar])*(1+Rhostar)*(1 - muestar*Gstar*(1+zeta_gw/zeta_zw/Rhostar)) - m[:gammstar]/m[:bet]*(1+Rhostar)
    m[:zeta_nn]     = m[:gammstar]/m[:bet] + m[:gammstar]*Rkstar/m[:pistar]/exp(m[:zstar])*(1+Rhostar)*muestar*Gstar*zeta_gw/zeta_zw/Rhostar
    m[:zeta_nmue]   = m[:gammstar]*Rkstar/m[:pistar]/exp(m[:zstar])*(1+Rhostar)*muestar*Gstar*(1 - zeta_gw*zeta_zmue/zeta_zw)
    m[:zeta_nsigw]  = m[:gammstar]*Rkstar/m[:pistar]/exp(m[:zstar])*(1+Rhostar)*muestar*Gstar*(zeta_Gsigw-zeta_gw/zeta_zw*zeta_zsigw)

    return m
end
