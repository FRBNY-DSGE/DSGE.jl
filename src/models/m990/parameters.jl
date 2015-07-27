# Then assign parameters to a theta vector
# θ = Parameters(α, β, etc.)
type Parameters990 <: Parameters
    alp::Param               # α     
    zeta_p::Param            # ζ_p   
    iota_p::Param            # ι_p   
    del::Param               # δ     
    ups::Param               # υ     
    Bigphi::Param            # Φ     
    s2::Param                # s2    
    h::Param                 # h     
    ppsi::Param              # ppsi  
    nu_l::Param              # ν_l   
    zeta_w::Param            # ζ_w   
    iota_w::Param            # ι_w   
    law::Param               # λ_w   
    bet::Param               # β     
    psi1::Param              # ψ₁    
    psi2::Param              # ψ₂    
    psi3::Param              # ψ₃    
    pistar::Param            # π_star
    sigmac::Param            # σ_c   
    rho::Param               # ρ     
    epsp::Param              # ϵ_p   
    epsw::Param              # ϵ_w   
    Fom::Param               # Fω    
    sprd::Param              # sprd  
    zeta_spb::Param          # ζ_spb 
    gammstar::Param          # γ_star
    gam::Param               # γ     
    Lmean::Param             # Lmean 
    gstar::Param             # g_star
    ρ_g::Param               # ρ_g   
    ρ_b::Param               # ρ_b   
    ρ_mu::Param              # ρ_μ   
    ρ_z::Param               # ρ_z   
    ρ_laf::Param             # ρ_λ_f 
    ρ_law::Param             # ρ_λ_w 
    ρ_rm::Param              # ρ_rm  
    ρ_sigw::Param            # ρ_σ_w 
    ρ_mue::Param             # ρ_μ_e 
    ρ_gamm::Param            # ρ_γ   
    ρ_pist::Param            # ρ_π_star
    ρ_lr::Param              # ρ_lr    
    ρ_zp::Param              # ρ_zp    
    ρ_tfp::Param             # ρ_tfp   
    ρ_gdpdef::Param          # ρ_gdpdef
    ρ_pce::Param             # ρ_pce   
    σ_g::Param               # σ_g     
    σ_b::Param               # σ_b     
    σ_mu::Param              # σ_μ     
    σ_z::Param               # σ_z     
    σ_laf::Param             # σ_λ_f   
    σ_law::Param             # σ_λ_w   
    σ_rm::Param              # σ_rm    
    σ_sigw::Param            # σ_σ_w   
    σ_mue::Param             # σ_μ_e   
    σ_gamm::Param            # σ_γ     
    σ_pist::Param            # σ_π_star
    σ_lr::Param              # σ_lr    
    σ_zp::Param              # σ_zp    
    σ_tfp::Param             # σ_tfp   
    σ_gdpdef::Param          # σ_gdpdef
    σ_pce::Param             # σ_pce   
    σ_rm1::Param             # σ_rm1   
    σ_rm2::Param             # σ_rm2   
    σ_rm3::Param             # σ_rm3   
    σ_rm4::Param             # σ_rm4   
    σ_rm5::Param             # σ_rm5   
    σ_rm6::Param             # σ_rm6   
    σ_rm7::Param             # σ_rm7   
    σ_rm8::Param             # σ_rm8   
    σ_rm9::Param             # σ_rm9   
    σ_rm10::Param            # σ_rm10  
    σ_rm11::Param            # σ_rm11  
    σ_rm12::Param            # σ_rm12  
    σ_rm13::Param            # σ_rm13  
    σ_rm14::Param            # σ_rm14  
    σ_rm15::Param            # σ_rm15  
    σ_rm16::Param            # σ_rm16  
    σ_rm17::Param            # σ_rm17  
    σ_rm18::Param            # σ_rm18  
    σ_rm19::Param            # σ_rm19  
    σ_rm20::Param            # σ_rm20  
    eta_gz::Param            # η_gz    
    eta_laf::Param           # η_λ_f   
    eta_law::Param           # η_λ_w   
    modelalp_ind::Param      # modelα_ind
    gamm_gdpdef::Param       # γ_gdpdef  
    del_gdpdef::Param        # δ_gdpdef  

    zstar::Float64           # z_star  
    rstar::Float64           # r_star 
    Rstarn::Float64          # R_starn
    rkstar::Float64          # rk_star
    wstar::Float64           # w_star 
    Lstar::Float64           # L_star 
    kstar::Float64           # k_star 
    kbarstar::Float64        # k̄_star 
    istar::Float64           # i_star 
    ystar::Float64           # y_star 
    cstar::Float64           # c_star 
    wl_c::Float64            # wl_c   
    nstar::Float64           # n_star  
    vstar::Float64           # v_star  
    zeta_spsigw::Float64     # ζ_spsigw
    zeta_spmue::Float64      # ζ_spmue 
    zeta_nRk::Float64        # ζ_nRk   
    zeta_nR::Float64         # ζ_nR    
    zeta_nqk::Float64        # ζ_nqk   
    zeta_nn::Float64         # ζ_nn    
    zeta_nmue::Float64       # ζ_nmue  
    zeta_nsigw::Float64      # ζ_nsigw 
    
    # This constructor takes in as arguments and initializes the Param fields, then
    #   immediately calls steadystate!() to calculate and initialize the steady-state values
    #   (which have type Float64)
    function Parameters990(alp, zeta_p, iota_p, del, ups, Bigphi, s2, h, ppsi, nu_l, zeta_w,
    iota_w, law, bet, psi1, psi2, psi3, pistar, sigmac, rho, epsp, epsw, Fom, sprd,
    zeta_spb, gammstar, gam, Lmean, gstar, ρ_g, ρ_b, ρ_mu, ρ_z, ρ_laf, ρ_law, ρ_rm,
    ρ_sigw, ρ_mue, ρ_gamm, ρ_pistar, ρ_lr, ρ_zp, ρ_tfp, ρ_gdpdef, ρ_pce, σ_g, σ_b, σ_mu,
    σ_z, σ_laf, σ_law, σ_rm, σ_sigw, σ_mue, σ_gamm, σ_pistar, σ_lr, σ_zp, σ_tfp, σ_gdpdef,
    σ_pce, σ_rm1, σ_rm2, σ_rm3, σ_rm4, σ_rm5, σ_rm6, σ_rm7, σ_rm8, σ_rm9, σ_rm10, σ_rm11,
    σ_rm12, σ_rm13, σ_rm14, σ_rm15, σ_rm16, σ_rm17, σ_rm18, σ_rm19, σ_rm20, eta_gz,
    eta_laf, eta_law, modelalp_ind, gamm_gdpdef, del_gdpdef)
        steadystate!(new(alp, zeta_p, iota_p, del, ups, Bigphi, s2, h, ppsi, nu_l, zeta_w,
        iota_w, law, bet, psi1, psi2, psi3, pistar, sigmac, rho, epsp, epsw, Fom, sprd,
        zeta_spb, gammstar, gam, Lmean, gstar, ρ_g, ρ_b, ρ_mu, ρ_z, ρ_laf, ρ_law, ρ_rm,
        ρ_sigw, ρ_mue, ρ_gamm, ρ_pistar, ρ_lr, ρ_zp, ρ_tfp, ρ_gdpdef, ρ_pce, σ_g, σ_b, σ_mu,
        σ_z, σ_laf, σ_law, σ_rm, σ_sigw, σ_mue, σ_gamm, σ_pistar, σ_lr, σ_zp, σ_tfp, σ_gdpdef,
        σ_pce, σ_rm1, σ_rm2, σ_rm3, σ_rm4, σ_rm5, σ_rm6, σ_rm7, σ_rm8, σ_rm9, σ_rm10, σ_rm11,
        σ_rm12, σ_rm13, σ_rm14, σ_rm15, σ_rm16, σ_rm17, σ_rm18, σ_rm19, σ_rm20, eta_gz,
        eta_laf, eta_law, modelalp_ind, gamm_gdpdef, del_gdpdef))
    end
end



# TODO: some parameters (e.g. s2) have type = 0 but a and b
# Instantiate Parameters990 type
function Parameters990(spec::Dict{String, Any})
    alp = Param(0.1596, false, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999))
    zeta_p = Param(0.8940, false, (1e-5, 0.999), Beta(0.5, 0.1), 1, (1e-5, 0.999))
    iota_p = Param(0.1865, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999))
    del = Param(0.025)
    ups = Param(1.000, true, (0., 10.), Gamma(1., 0.5), 2, (1e-5, 0.))
    Bigphi = Param(1.1066, false, (1., 10.), Normal(1.25, 0.12), 2, (1.00, 10.00))
    s2 = Param(2.7314, false, (-15., 15.), Normal(4., 1.5), 0, (-15., 15.))
    h = Param(0.5347, false, (1e-5, 0.999), Beta(0.7, 0.1), 1, (1e-5, 0.999))
    ppsi = Param(0.6862, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999))
    nu_l = Param(2.5975, false, (1e-5, 10.), Normal(2, 0.75), 2, (1e-5, 10.))
    zeta_w = Param(0.9291, false, (1e-5, 0.999), Beta(0.5, 0.1), 1, (1e-5, 0.999))
    iota_w = Param(0.2992, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999))
    law = Param(1.5)
    # laf = []
    bet = Param(0.1402, scalefunction = x -> 1/(1 + x/100), false, (1e-5, 10.), Gamma(0.25, 0.1), 2, (1e-5, 10.))
    psi1 = Param(1.3679, false, (1e-5, 10.), Normal(1.5, 0.25), 2, (1e-5, 10.00))
    psi2 = Param(0.0388, false, (-0.5, 0.5), Normal(0.12, 0.05), 0, (-0.5, 0.5))
    psi3 = Param(0.2464, false, (-0.5, 0.5), Normal(0.12, 0.05), 0, (-0.5, 0.5))
    pistar = Param(0.5000, scalefunction = x -> 1 + x/100, true, (1e-5, 10.), Gamma(0.75, 0.4), 2, (1e-5, 10.))
    sigmac = Param(0.8719, false, (1e-5, 10.), Normal(1.5, 0.37), 2, (1e-5, 10.))
    rho = Param(0.7126, false, (1e-5, 0.999), Beta(0.75, 0.10), 1, (1e-5, 0.999))
    epsp = Param(10.)
    epsw = Param(10.)

    # financial frictions parameters
    Fom = Param(0.0300, scalefunction = x -> 1 - (1-x)^0.25, true, (1e-5, 0.99999), Beta(0.03, 0.01), 1, (1e-5, 0.99))
    sprd = Param(1.7444, scalefunction = x -> (1 + x/100)^0.25, false, (0., 100.), Gamma(2., 0.1), 2, (1e-5, 0.))
    zeta_spb = Param(0.0559, false, (1e-5, 0.99999), Beta(0.05, 0.005), 1, (1e-5, 0.99))
    gammstar = Param(0.9900, true, (1e-5, 0.99999), Beta(0.99, 0.002), 1, (1e-5, 0.99))

    # exogenous processes - level
    gam = Param(0.3673, scalefunction = x -> x/100, false, (-5., 5.), Normal(0.4, 0.1), 0, (-5.0, 5.0))
    Lmean = Param(-45.9364, false, (-1000., 1000.), Normal(-45, 5), 0, (-1000., 1000.))
    gstar = Param(0.18)

    # exogenous processes - autocorrelation
    ρ_g = Param(0.9863, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_b = Param(0.9410, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_mu = Param(0.8735, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_z = Param(0.9446, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_laf = Param(0.8827, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_law = Param(0.3884, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_rm = Param(0.2135, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_sigw = Param(0.9898, false, (1e-5, 0.99999), Beta(0.75, 0.15), 1, (1e-5, 0.99))
    ρ_mue = Param(0.7500, true, (1e-5, 0.99999), Beta(0.75, 0.15), 1, (1e-5, 0.99))
    ρ_gamm = Param(0.7500, true, (1e-5, 0.99999), Beta(0.75, 0.15), 1, (1e-5, 0.99))
    ρ_pist = Param(0.9900, true, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_lr = Param(0.6936, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_zp = Param(0.8910, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_tfp = Param(0.1953, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_gdpdef = Param(0.5379, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))
    ρ_pce = Param(0.2320, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999))

    # exogenous processes - standard deviation
    σ_g = Param(2.5230, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_b = Param(0.0292, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_mu = Param(0.4559, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_z = Param(0.6742, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_laf = Param(0.1314, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_law = Param(0.3864, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_rm = Param(0.2380, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_sigw = Param(0.0428, false, (1e-7, 100.), RootInverseGamma(4.00, 0.05), 2, (1e-5, 0.))
    σ_mue = Param(0., true, (1e-7, 100.), RootInverseGamma(4.00, 0.05), 2, (1e-5, 0.))
    σ_gamm = Param(0., true, (1e-7, 100.), RootInverseGamma(4.00, 0.01), 2, (1e-5, 0.))
    σ_pist = Param(0.0269, false, (1e-8, 5.), RootInverseGamma(6., 0.03), 2, (1e-8, 5.))
    σ_lr = Param(0.1766, false, (1e-8, 10.), RootInverseGamma(2., 0.75), 2, (1e-8, 5.))
    σ_zp = Param(0.1662, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_tfp = Param(0.9391, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_gdpdef = Param(0.1575, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))
    σ_pce = Param(0.0999, false, (1e-8, 5.), RootInverseGamma(2.00, 0.1), 2, (1e-8, 5.))

    # standard deviations of the anticipated policy shocks
    for i = 1:spec["nantpad"]
        if i < 13
            eval(parse("σ_rm$i = Param(0.20, false, (1e-7, 100.), RootInverseGamma(4.00, 0.2), 2, (1e-5, 0.))"))
        else
            eval(parse("σ_rm$i = Param(0.0, true, (1e-7, 100.), RootInverseGamma(4.00, 0.2), 2, (1e-5, 0.))"))
        end
    end

    eta_gz = Param(0.8400, false, (1e-5, 0.999), Beta(0.50, 0.20), 1, (1e-5, 0.999))
    eta_laf = Param(0.7892, false, (1e-5, 0.999), Beta(0.50, 0.20), 1, (1e-5, 0.999))
    eta_law = Param(0.4226, false, (1e-5, 0.999), Beta(0.50, 0.20), 1, (1e-5, 0.999))

    modelalp_ind = Param(0.0000, true, (0.000, 1.000), Beta(0.50, 0.20), 0, (0., 0.))
    gamm_gdpdef = Param(1.0354, false, (-10., 10.), Normal(1.00, 2.), 0, (-10., -10.))
    del_gdpdef = Param(0.0181, false, (-9.1, 9.1), Normal(0.00, 2.), 0, (-10., -10.))


    
    Θ = Parameters990(alp, zeta_p, iota_p, del, ups, Bigphi, s2, h, ppsi, nu_l, zeta_w,
                      iota_w, law, bet, psi1, psi2, psi3, pistar, sigmac, rho, epsp,
                      epsw, Fom, sprd, zeta_spb, gammstar, gam, Lmean, gstar, ρ_g, ρ_b,
                      ρ_mu, ρ_z, ρ_laf, ρ_law, ρ_rm, ρ_sigw, ρ_mue, ρ_gamm, ρ_pist, ρ_lr,
                      ρ_zp, ρ_tfp, ρ_gdpdef, ρ_pce, σ_g, σ_b, σ_mu, σ_z, σ_laf, σ_law,
                      σ_rm, σ_sigw, σ_mue, σ_gamm, σ_pist, σ_lr, σ_zp, σ_tfp, σ_gdpdef,
                      σ_pce, σ_rm1, σ_rm2, σ_rm3, σ_rm4, σ_rm5, σ_rm6, σ_rm7, σ_rm8,
                      σ_rm9, σ_rm10, σ_rm11, σ_rm12, σ_rm13, σ_rm14, σ_rm15, σ_rm16,
                      σ_rm17, σ_rm18, σ_rm19, σ_rm20, eta_gz, eta_laf, eta_law,
                      modelalp_ind, gamm_gdpdef, del_gdpdef)
    spec["n_params"] = count(field -> isa(getfield(Θ, field), Param), names(Θ))
    return Θ
end



# (Re)calculates steady-state values from Param fields in Θ
# The functions to calculate financial frictions additions (ζ_spb_fn, etc.) can be found in
# init/FinancialFrictionsFunctions.jl
function steadystate!(Θ::Parameters990)
    Θ.zstar = log(1+Θ.gam) + Θ.alp/(1-Θ.alp)*log(Θ.ups)
    Θ.rstar = exp(Θ.sigmac*Θ.zstar) / Θ.bet
    Θ.Rstarn = 100*(Θ.rstar*Θ.pistar - 1)
    Θ.rkstar = Θ.sprd*Θ.rstar*Θ.ups - (1-Θ.del)
    Θ.wstar = (Θ.alp^Θ.alp * (1-Θ.alp)^(1-Θ.alp) * Θ.rkstar^(-Θ.alp) / Θ.Bigphi)^(1/(1-Θ.alp))
    Θ.Lstar = 1.
    Θ.kstar = (Θ.alp/(1-Θ.alp)) * Θ.wstar * Θ.Lstar / Θ.rkstar
    Θ.kbarstar = Θ.kstar * (1+Θ.gam) * Θ.ups^(1 / (1-Θ.alp))
    Θ.istar = Θ.kbarstar * (1-((1-Θ.del)/((1+Θ.gam) * Θ.ups^(1/(1-Θ.alp)))))
    Θ.ystar = (Θ.kstar^Θ.alp) * (Θ.Lstar^(1-Θ.alp)) / Θ.Bigphi
    Θ.cstar = (1-Θ.gstar)*Θ.ystar - Θ.istar
    Θ.wl_c = (Θ.wstar*Θ.Lstar)/(Θ.cstar*Θ.law)

    # FINANCIAL FRICTIONS ADDITIONS
    # solve for sigmaomegastar and zomegastar
    zwstar = quantile(Normal(), Θ.Fom.scaledvalue)
    sigwstar = fzero(sigma -> ζ_spb_fn(zwstar, sigma, Θ.sprd) - Θ.zeta_spb, 0.5)

    # evaluate omegabarstar
    omegabarstar = ω_fn(zwstar, sigwstar)

    # evaluate all BGG function elasticities
    Gstar = G_fn(zwstar, sigwstar)
    Gammastar = Γ_fn(zwstar, sigwstar)
    dGdomegastar = dG_dω_fn(zwstar, sigwstar)
    d2Gdomega2star = d2G_dω2_fn(zwstar, sigwstar)
    dGammadomegastar = dΓ_dω_fn(zwstar)
    d2Gammadomega2star = d2Γ_dω2_fn(zwstar, sigwstar)
    dGdsigmastar = dG_dσ_fn(zwstar, sigwstar)
    d2Gdomegadsigmastar = d2G_dωdσ_fn(zwstar, sigwstar)
    dGammadsigmastar = dΓ_dσ_fn(zwstar, sigwstar)
    d2Gammadomegadsigmastar = d2Γ_dωdσ_fn(zwstar, sigwstar)

    # evaluate mu, nk, and Rhostar
    muestar = μ_fn(zwstar, sigwstar, Θ.sprd)
    nkstar = nk_fn(zwstar, sigwstar, Θ.sprd)
    Rhostar = 1/nkstar - 1

    # evaluate wekstar and vkstar
    wekstar = (1-Θ.gammstar/Θ.bet)*nkstar - Θ.gammstar/Θ.bet*(Θ.sprd*(1-muestar*Gstar) - 1)
    vkstar = (nkstar-wekstar)/Θ.gammstar

    # evaluate nstar and vstar
    Θ.nstar = nkstar*Θ.kstar
    Θ.vstar = vkstar*Θ.kstar

    # a couple of combinations
    GammamuG = Gammastar - muestar*Gstar
    GammamuGprime = dGammadomegastar - muestar*dGdomegastar

    # elasticities wrt omegabar
    zeta_bw = ζ_bω_fn(zwstar, sigwstar, Θ.sprd)
    zeta_zw = ζ_zω_fn(zwstar, sigwstar, Θ.sprd)
    zeta_bw_zw = zeta_bw/zeta_zw

    # elasticities wrt sigw
    zeta_bsigw = sigwstar * (((1 - muestar*dGdsigmastar/dGammadsigmastar) / (1 - muestar*dGdomegastar/dGammadomegastar) - 1)*dGammadsigmastar*Θ.sprd + muestar*nkstar*(dGdomegastar*d2Gammadomegadsigmastar - dGammadomegastar*d2Gdomegadsigmastar)/GammamuGprime^2) / ((1 - Gammastar)*Θ.sprd + dGammadomegastar/GammamuGprime*(1-nkstar))
    zeta_zsigw = sigwstar * (dGammadsigmastar - muestar*dGdsigmastar) / GammamuG
    Θ.zeta_spsigw = (zeta_bw_zw*zeta_zsigw - zeta_bsigw) / (1-zeta_bw_zw)
    
    # elasticities wrt mue
    zeta_bmue = muestar * (nkstar*dGammadomegastar*dGdomegastar/GammamuGprime+dGammadomegastar*Gstar*Θ.sprd) / ((1-Gammastar)*GammamuGprime*Θ.sprd + dGammadomegastar*(1-nkstar))
    zeta_zmue = -muestar*Gstar/GammamuG
    Θ.zeta_spmue = (zeta_bw_zw*zeta_zmue - zeta_bmue) / (1-zeta_bw_zw)

    # some ratios/elasticities
    Rkstar = Θ.sprd*Θ.pistar*Θ.rstar # (rkstar+1-delta)/ups*pistar
    zeta_gw = dGdomegastar/Gstar*omegabarstar
    zeta_Gsigw = dGdsigmastar/Gstar*sigwstar
    
    # elasticities for the net worth evolution
    Θ.zeta_nRk = Θ.gammstar*Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Rhostar)*(1 - muestar*Gstar*(1 - zeta_gw/zeta_zw))
    Θ.zeta_nR = Θ.gammstar/Θ.bet*(1+Rhostar)*(1 - nkstar + muestar*Gstar*Θ.sprd*zeta_gw/zeta_zw)
    Θ.zeta_nqk = Θ.gammstar*Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Rhostar)*(1 - muestar*Gstar*(1+zeta_gw/zeta_zw/Rhostar)) - Θ.gammstar/Θ.bet*(1+Rhostar)
    Θ.zeta_nn = Θ.gammstar/Θ.bet + Θ.gammstar*Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Rhostar)*muestar*Gstar*zeta_gw/zeta_zw/Rhostar
    Θ.zeta_nmue = Θ.gammstar*Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Rhostar)*muestar*Gstar*(1 - zeta_gw*zeta_zmue/zeta_zw)
    Θ.zeta_nsigw = Θ.gammstar*Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Rhostar)*muestar*Gstar*(zeta_Gsigw-zeta_gw/zeta_zw*zeta_zsigw)
    
    return Θ
end
