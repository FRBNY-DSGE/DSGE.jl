using Distributions
using DistributionsExt: Beta, Gamma, InverseGamma

include("../../code/abstractmodel/param.jl")
include("../../code/abstractmodel/parameters.jl")
include("../../code/init/ff_fcns.jl")

# Then assign parameters to a theta vector
# θ = Parameters(α, β, etc.)
type Parameters990 <: Parameters
    alp::Param
    zeta_p::Param
    iota_p::Param
    del::Param
    ups::Param
    Bigphi::Param
    s2::Param
    h::Param
    ppsi::Param
    nu_l::Param
    zeta_w::Param
    iota_w::Param
    law::Param
    bet::Param
    psi1::Param
    psi2::Param
    psi3::Param
    pistar::Param
    sigmac::Param
    rho::Param
    epsp::Param
    epsw::Param
    Fom::Param
    sprd::Param
    zeta_spb::Param
    gammstar::Param
    gam::Param
    Lmean::Param
    gstar::Param
    ρ_g::Param
    ρ_b::Param
    ρ_mu::Param
    ρ_z::Param
    ρ_laf::Param
    ρ_law::Param
    ρ_rm::Param
    ρ_sigw::Param
    ρ_mue::Param
    ρ_gamm::Param
    ρ_pist::Param
    ρ_lr::Param
    ρ_zp::Param
    ρ_tfp::Param
    ρ_gdpdef::Param
    ρ_pce::Param
    σ_g::Param
    σ_b::Param
    σ_mu::Param
    σ_z::Param
    σ_laf::Param
    σ_law::Param
    σ_rm::Param
    σ_sigw::Param
    σ_mue::Param
    σ_gamm::Param
    σ_pist::Param
    σ_lr::Param
    σ_zp::Param
    σ_tfp::Param
    σ_gdpdef::Param
    σ_pce::Param
    σ_rm1::Param
    σ_rm2::Param
    σ_rm3::Param
    σ_rm4::Param
    σ_rm5::Param
    σ_rm6::Param
    σ_rm7::Param
    σ_rm8::Param
    σ_rm9::Param
    σ_rm10::Param
    σ_rm11::Param
    σ_rm12::Param
    σ_rm13::Param
    σ_rm14::Param
    σ_rm15::Param
    σ_rm16::Param
    σ_rm17::Param
    σ_rm18::Param
    σ_rm19::Param
    σ_rm20::Param
    eta_gz::Param
    eta_laf::Param
    eta_law::Param
    modelalp_ind::Param
    gamm_gdpdef::Param
    del_gdpdef::Param

    zstar::Float64
    rstar::Float64
    Rstarn::Float64
    rkstar::Float64
    wstar::Float64
    Lstar::Float64
    kstar::Float64
    kbarstar::Float64
    istar::Float64
    ystar::Float64
    cstar::Float64
    wl_c::Float64
    zwstar::Float64
    sigwstar::Float64
    omegabarstar::Float64
    Gstar::Float64
    Gammastar::Float64
    dGdomegastar::Float64
    d2Gdomega2star::Float64
    dGammadomegastar::Float64
    d2Gammadomega2star::Float64
    dGdsigmastar::Float64
    d2Gdomegadsigmastar::Float64
    dGammadsigmastar::Float64
    d2Gammadomegadsigmastar::Float64
    muestar::Float64
    nkstar::Float64
    Rhostar::Float64
    wekstar::Float64
    vkstar::Float64
    nstar::Float64
    vstar::Float64
    GammamuG::Float64
    GammamuGprime::Float64
    zeta_bw::Float64
    zeta_zw::Float64
    zeta_bw_zw::Float64
    zeta_bsigw::Float64
    zeta_zsigw::Float64
    zeta_spsigw::Float64
    zeta_bmue::Float64
    zeta_zmue::Float64
    zeta_spmue::Float64
    Rkstar::Float64
    zeta_Gw::Float64
    zeta_Gsigw::Float64
    zeta_nRk::Float64
    zeta_nR::Float64
    zeta_nqk::Float64
    zeta_nn::Float64
    zeta_nmue::Float64
    zeta_nsigw::Float64
    
    # This constructor takes in as arguments and initializes the Param fields, then immediately calls steadystate() to
    # calculate and initialize the steady-state values (which have type Float64)
    function Parameters990(alp, zeta_p, iota_p, del, ups, Bigphi, s2, h, ppsi, nu_l, zeta_w, iota_w, law, bet, psi1, psi2, psi3,
    pistar, sigmac, rho, epsp, epsw, Fom, sprd, zeta_spb, gammstar, gam, Lmean, gstar, ρ_g, ρ_b, ρ_mu, ρ_z, ρ_laf,
    ρ_law, ρ_rm, ρ_sigw, ρ_mue, ρ_gamm, ρ_pistar, ρ_lr, ρ_zp, ρ_tfp, ρ_gdpdef, ρ_pce, σ_g, σ_b, σ_mu, σ_z, σ_laf, σ_law,
    σ_rm, σ_sigw, σ_mue, σ_gamm, σ_pistar, σ_lr, σ_zp, σ_tfp, σ_gdpdef, σ_pce, σ_rm1, σ_rm2, σ_rm3, σ_rm4, σ_rm5, σ_rm6,
    σ_rm7, σ_rm8, σ_rm9, σ_rm10, σ_rm11, σ_rm12, σ_rm13, σ_rm14, σ_rm15, σ_rm16, σ_rm17, σ_rm18, σ_rm19, σ_rm20, eta_gz,
    eta_laf, eta_law, modelalp_ind, gamm_gdpdef, del_gdpdef) 
      steadystate(new(alp, zeta_p, iota_p, del, ups, Bigphi,
      s2, h, ppsi, nu_l, zeta_w, iota_w, law, bet, psi1, psi2, psi3, pistar, sigmac, rho, epsp, epsw, Fom, sprd, zeta_spb,
      gammstar, gam, Lmean, gstar, ρ_g, ρ_b, ρ_mu, ρ_z, ρ_laf, ρ_law, ρ_rm, ρ_sigw, ρ_mue, ρ_gamm, ρ_pistar, ρ_lr, ρ_zp,
      ρ_tfp, ρ_gdpdef, ρ_pce, σ_g, σ_b, σ_mu, σ_z, σ_laf, σ_law, σ_rm, σ_sigw, σ_mue, σ_gamm, σ_pistar, σ_lr, σ_zp, σ_tfp,
      σ_gdpdef, σ_pce, σ_rm1, σ_rm2, σ_rm3, σ_rm4, σ_rm5, σ_rm6, σ_rm7, σ_rm8, σ_rm9, σ_rm10, σ_rm11, σ_rm12, σ_rm13,
      σ_rm14, σ_rm15, σ_rm16, σ_rm17, σ_rm18, σ_rm19, σ_rm20, eta_gz, eta_laf, eta_law, modelalp_ind, gamm_gdpdef,
      del_gdpdef))
    end

end



# TODO: some parameters (e.g. s2) have type = 0 but a and b
# Instantiate Parameters990 type
function Parameters990()
    alp = Param(0.1596, false, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999)) # alp
    zeta_p = Param(0.8940, false, (1e-5, 0.999), Beta(0.5, 0.1), 1, (1e-5, 0.999)) # zeta_p
    iota_p = Param(0.1865, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999)) # iota_p
    del = Param(0.025) # del = 0.025
    ups = Param(1.000, true, (0., 10.), Gamma(1., 0.5), 2, (1e-5, 0.)) # ups
    Bigphi = Param(1.1066, false, (1., 10.), Normal(1.25, 0.12), 2, (1.00, 10.00)) # Bigphi
    s2 = Param(2.7314, false, (-15., 15.), Normal(4., 1.5), 0, (-15., 15.)) # s2
    h = Param(0.5347, false, (1e-5, 0.999), Beta(0.7, 0.1), 1, (1e-5, 0.999)) # h
    ppsi = Param(0.6862, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999)) # ppsi
    nu_l = Param(2.5975, false, (1e-5, 10.), Normal(2, 0.75), 2, (1e-5, 10.)) # nu_l
    zeta_w = Param(0.9291, false, (1e-5, 0.999), Beta(0.5, 0.1), 1, (1e-5, 0.999)) # zeta_w
    iota_w = Param(0.2992, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999)) # iota_w
    law = Param(1.5) # law = 1.5;
    # laf = [];
    bet = Param(0.1402, transf = x -> 1/(1 + x/100), false, (1e-5, 10.), Gamma(0.25, 0.1), 2, (1e-5, 10.)) # bet
    psi1 = Param(1.3679, false, (1e-5, 10.), Normal(1.5, 0.25), 2, (1e-5, 10.00)) # psi1
    psi2 = Param(0.0388, false, (-0.5, 0.5), Normal(0.12, 0.05), 0, (-0.5, 0.5)) # psi2
    psi3 = Param(0.2464, false, (-0.5, 0.5), Normal(0.12, 0.05), 0, (-0.5, 0.5)) # psi3
    pistar = Param(0.5000, transf = x -> 1 + x/100, false, (1e-5, 10.), Gamma(0.75, 0.4), 2, (1e-5, 10.)) # pistar
    sigmac = Param(0.8719, false, (1e-5, 10.), Normal(1.5, 0.37), 2, (1e-5, 10.)) # sigmac
    rho = Param(0.7126, false, (1e-5, 0.999), Beta(0.75, 0.10), 1, (1e-5, 0.999)) # rho
    epsp = Param(10.) # epsp
    epsw = Param(10.) # epsw

    nantpad = 20
    
    # financial frictions parameters
    Fom = Param(0.0300, transf = x -> 1 - (1-x)^0.25, true, (1e-5, 0.999), Beta(0.03, 0.01), 1, (1e-5, 0.999)) # Fom
    sprd = Param(1.7444, transf = x -> (1 + x/100)^0.25, false, (0., 100.), Gamma(2., 0.1), 2, (1e-5, 0.)) # sprd
    zeta_spb = Param(0.0559, false, (1e-5, 0.999), Beta(0.05, 0.005), 1, (1e-5, 0.999)) # zeta_spb
    gammstar = Param(0.9900, true, (1e-5, 0.999), Beta(0.99, 0.002), 1, (1e-5, 0.999))  # gammstar

    # exogenous processes - level
    gam = Param(0.3673, transf = x -> x/100, false, (-5., 5.), Normal(0.4, 0.1), 0, (-5.0, 5.0)) # gam
    Lmean = Param(-45.9364, false, (-1000., 1000.), Normal(-45, 5), 0, (-1000., 1000.)) # Lmean
    gstar = Param(0.18) # gstar = .18;

    # exogenous processes - autocorrelation
    ρ_g = Param(0.9863, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_g
    ρ_b = Param(0.9410, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_b
    ρ_mu = Param(0.8735, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_mu
    ρ_z = Param(0.9446, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_z
    ρ_laf = Param(0.8827, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_laf
    ρ_law = Param(0.3884, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_law
    ρ_rm = Param(0.2135, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_rm
    ρ_sigw = Param(0.9898, false, (1e-5, 0.999), Beta(0.75, 0.15), 1, (1e-5, 0.999)) # rho_sigw
    ρ_mue = Param(0.7500, true, (1e-5, 0.999), Beta(0.75, 0.15), 1, (1e-5, 0.999)) # rho_mue
    ρ_gamm = Param(0.7500, true, (1e-5, 0.999), Beta(0.75, 0.15), 1, (1e-5, 0.999)) # rho_gamm
    ρ_pist = Param(0.9900, true, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_pist
    ρ_lr = Param(0.6936, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_lr
    ρ_zp = Param(0.8910, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_zp
    ρ_tfp = Param(0.1953, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_tfp
    ρ_gdpdef = Param(0.5379, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_gdpdef
    ρ_pce = Param(0.2320, false, (1e-5, 0.999), Beta(0.5, 0.2), 1, (1e-5, 0.999)) # rho_pce

    # exogenous processes - standard deviation
    σ_g = Param(2.5230, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_g
    σ_b = Param(0.0292, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_b
    σ_mu = Param(0.4559, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_mu
    σ_z = Param(0.6742, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_z
    σ_laf = Param(0.1314, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_laf
    σ_law = Param(0.3864, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_law
    σ_rm = Param(0.2380, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_rm
    σ_sigw = Param(0.0428, false, (1e-7, 100.), InverseGamma(0.05, 4.00), 2, (1e-5, 0.)) # sig_sigw
    σ_mue = Param(0., true, (1e-7, 100.), InverseGamma(0.05, 4.00), 2, (1e-5, 0.)) # sig_mue
    σ_gamm = Param(0., true, (1e-7, 100.), InverseGamma(0.01, 4.00), 2, (1e-5, 0.)) # sig_gamm
    σ_pist = Param(0.0269, false, (1e-8, 5.), InverseGamma(0.03, 6.), 2, (1e-8, 5.)) # sig_pist
    σ_lr = Param(0.1766, false, (1e-8, 10.), InverseGamma(0.75, 2.), 2, (1e-8, 5.)) # sig_lr
    σ_zp = Param(0.1662, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_zp
    σ_tfp = Param(0.9391, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_tfp
    σ_gdpdef = Param(0.1575, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_gdpdef
    σ_pce = Param(0.0999, false, (1e-8, 5.), InverseGamma(0.1, 2.00), 2, (1e-8, 5.)) # sig_pce

    # TODO: better way of representing these?
    # standard deviations of the anticipated policy shocks
    for i = 1:nantpad
        if i < 13
            eval(parse("σ_rm$i = Param(0.20, false, (1e-7, 100.), InverseGamma(0.2, 4.00), 2, (1e-5, 0.))"))
        else
            eval(parse("σ_rm$i = Param(0.20, true, (1e-7, 100.), InverseGamma(0.2, 4.00), 2, (1e-5, 0.))"))
        end
    end

    eta_gz = Param(0.8400, false, (1e-5, 0.999), Beta(0.50, 0.20), 1, (1e-5, 0.999)) # eta_gz
    eta_laf = Param(0.7892, false, (1e-5, 0.999), Beta(0.50, 0.20), 1, (1e-5, 0.999)) # eta_laf
    eta_law = Param(0.4226, false, (1e-5, 0.999), Beta(0.50, 0.20), 1, (1e-5, 0.999)) # eta_law

    modelalp_ind = Param(0.0000, true, (0.000, 1.000), Beta(0.50, 0.20), 0, (0., 0.)) # modelalp_ind
    gamm_gdpdef = Param(1.0354, false, (-10., 10.), Normal(1.00, 2.), 0, (-10., -10.)) # gamm_gdpdef
    del_gdpdef = Param(0.0181, false, (-9.1, 9.1), Normal(0.00, 2.), 0, (-10., -10.)) # del_gdpdef
    
    return Parameters990(alp, zeta_p, iota_p, del, ups, Bigphi, s2, h, ppsi, nu_l, zeta_w, iota_w, law, bet, psi1, psi2, psi3, pistar, sigmac, rho, epsp, epsw, Fom, sprd, zeta_spb, gammstar, gam, Lmean, gstar, ρ_g, ρ_b, ρ_mu, ρ_z, ρ_laf, ρ_law, ρ_rm, ρ_sigw, ρ_mue, ρ_gamm, ρ_pist, ρ_lr, ρ_zp, ρ_tfp, ρ_gdpdef, ρ_pce, σ_g, σ_b, σ_mu, σ_z, σ_laf, σ_law, σ_rm, σ_sigw, σ_mue, σ_gamm, σ_pist, σ_lr, σ_zp, σ_tfp, σ_gdpdef, σ_pce, σ_rm1, σ_rm2, σ_rm3, σ_rm4, σ_rm5, σ_rm6, σ_rm7, σ_rm8, σ_rm9, σ_rm10, σ_rm11, σ_rm12, σ_rm13, σ_rm14, σ_rm15, σ_rm16, σ_rm17, σ_rm18, σ_rm19, σ_rm20, eta_gz, eta_laf, eta_law, modelalp_ind, gamm_gdpdef, del_gdpdef)
    #return Parameters990(alp, zeta_p, iota_p, ups, Bigphi, s2, h, ppsi, nu_l, zeta_w, iota_w, bet, psi1, psi2, psi3, pistar, sigmac, rho, Fom, sprd, zeta_spb, gammstar, gam, Lmean, ρ_g, ρ_b, ρ_mu, ρ_z, ρ_laf, ρ_law, ρ_rm, ρ_sigw, ρ_mue, ρ_gamm, ρ_pist, ρ_lr, ρ_zp, ρ_tfp, ρ_gdpdef, ρ_pce, σ_g, σ_b, σ_mu, σ_z, σ_laf, σ_law, σ_rm, σ_sigw, σ_mue, σ_gamm, σ_pist, σ_lr, σ_zp, σ_tfp, σ_gdpdef, σ_pce, σ_rm1, σ_rm2, σ_rm3, σ_rm4, σ_rm5, σ_rm6, σ_rm7, σ_rm8, σ_rm9, σ_rm10, σ_rm11, σ_rm12, σ_rm13, σ_rm14, σ_rm15, σ_rm16, σ_rm17, σ_rm18, σ_rm19, σ_rm20, eta_gz, eta_laf, eta_law, modelalp_ind, gamm_gdpdef, del_gdpdef)
end



# (Re)calculates steady-state values from Param fields in Θ
# The functions called to calculate financial frictions additions (zetaspbfcn, etc.) can be found in
# code/init/ff_fcns.jl
function steadystate(Θ::Parameters990)
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
    Θ.zwstar = quantile(Normal(), Θ.Fom.tval)
    Θ.sigwstar = Roots.fzero(sigma -> zetaspbfcn(Θ.zwstar, sigma, Θ.sprd) - Θ.zeta_spb, 0.5)

    # evaluate omegabarstar
    Θ.omegabarstar = omegafcn(Θ.zwstar, Θ.sigwstar)

    # evaluate all BGG function elasticities
    Θ.Gstar = Gfcn(Θ.zwstar, Θ.sigwstar)
    Θ.Gammastar = Gammafcn(Θ.zwstar, Θ.sigwstar)
    Θ.dGdomegastar = dGdomegafcn(Θ.zwstar, Θ.sigwstar)
    Θ.d2Gdomega2star = d2Gdomega2fcn(Θ.zwstar, Θ.sigwstar)
    Θ.dGammadomegastar = dGammadomegafcn(Θ.zwstar)
    Θ.d2Gammadomega2star = d2Gammadomega2fcn(Θ.zwstar, Θ.sigwstar)
    Θ.dGdsigmastar = dGdsigmafcn(Θ.zwstar, Θ.sigwstar)
    Θ.d2Gdomegadsigmastar = d2Gdomegadsigmafcn(Θ.zwstar, Θ.sigwstar)
    Θ.dGammadsigmastar = dGammadsigmafcn(Θ.zwstar, Θ.sigwstar)
    Θ.d2Gammadomegadsigmastar = d2Gammadomegadsigmafcn(Θ.zwstar, Θ.sigwstar)

    # evaluate mu, nk, and Rhostar
    Θ.muestar = mufcn(Θ.zwstar, Θ.sigwstar, Θ.sprd)
    Θ.nkstar = nkfcn(Θ.zwstar, Θ.sigwstar, Θ.sprd)
    Θ.Rhostar = 1/Θ.nkstar - 1

    # evaluate wekstar and vkstar
    Θ.wekstar = (1-Θ.gammstar/Θ.bet)*Θ.nkstar - Θ.gammstar/Θ.bet*(Θ.sprd*(1-Θ.muestar*Θ.Gstar) - 1)
    Θ.vkstar = (Θ.nkstar-Θ.wekstar)/Θ.gammstar

    # evaluate nstar and vstar
    Θ.nstar = Θ.nkstar*Θ.kstar
    Θ.vstar = Θ.vkstar*Θ.kstar

    # a couple of combinations
    Θ.GammamuG = Θ.Gammastar - Θ.muestar*Θ.Gstar
    Θ.GammamuGprime = Θ.dGammadomegastar - Θ.muestar*Θ.dGdomegastar

    # elasticities wrt omegabar
    Θ.zeta_bw = zetabomegafcn(Θ.zwstar, Θ.sigwstar, Θ.sprd)
    Θ.zeta_zw = zetazomegafcn(Θ.zwstar, Θ.sigwstar, Θ.sprd)
    Θ.zeta_bw_zw = Θ.zeta_bw/Θ.zeta_zw

    # elasticities wrt sigw
    Θ.zeta_bsigw = Θ.sigwstar * (((1 - Θ.muestar*Θ.dGdsigmastar/Θ.dGammadsigmastar) / (1 - Θ.muestar*Θ.dGdomegastar/Θ.dGammadomegastar) - 1)*Θ.dGammadsigmastar*Θ.sprd + Θ.muestar*Θ.nkstar*(Θ.dGdomegastar*Θ.d2Gammadomegadsigmastar - Θ.dGammadomegastar*Θ.d2Gdomegadsigmastar)/Θ.GammamuGprime^2) / ((1 - Θ.Gammastar)*Θ.sprd + Θ.dGammadomegastar/Θ.GammamuGprime*(1-Θ.nkstar))
    Θ.zeta_zsigw = Θ.sigwstar * (Θ.dGammadsigmastar - Θ.muestar*Θ.dGdsigmastar) / Θ.GammamuG
    Θ.zeta_spsigw = (Θ.zeta_bw_zw*Θ.zeta_zsigw - Θ.zeta_bsigw) / (1-Θ.zeta_bw_zw)
    
    # elasticities wrt mue
    Θ.zeta_bmue = Θ.muestar * (Θ.nkstar*Θ.dGammadomegastar*Θ.dGdomegastar/Θ.GammamuGprime+Θ.dGammadomegastar*Θ.Gstar*Θ.sprd) / ((1-Θ.Gammastar)*Θ.GammamuGprime*Θ.sprd + Θ.dGammadomegastar*(1-Θ.nkstar))
    Θ.zeta_zmue = -Θ.muestar*Θ.Gstar/Θ.GammamuG
    Θ.zeta_spmue = (Θ.zeta_bw_zw*Θ.zeta_zmue - Θ.zeta_bmue) / (1-Θ.zeta_bw_zw)

    # some ratios/elasticities
    Θ.Rkstar = Θ.sprd*Θ.pistar*Θ.rstar # (rkstar+1-delta)/ups*pistar
    Θ.zeta_Gw = Θ.dGdomegastar/Θ.Gstar*Θ.omegabarstar
    Θ.zeta_Gsigw = Θ.dGdsigmastar/Θ.Gstar*Θ.sigwstar
    
    # elasticities for the net worth evolution
    Θ.zeta_nRk = Θ.gammstar*Θ.Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Θ.Rhostar)*(1 - Θ.muestar*Θ.Gstar*(1 - Θ.zeta_Gw/Θ.zeta_zw))
    Θ.zeta_nR = Θ.gammstar/Θ.bet*(1+Θ.Rhostar)*(1 - Θ.nkstar + Θ.muestar*Θ.Gstar*Θ.sprd*Θ.zeta_Gw/Θ.zeta_zw)
    Θ.zeta_nqk = Θ.gammstar*Θ.Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Θ.Rhostar)*(1 - Θ.muestar*Θ.Gstar*(1+Θ.zeta_Gw/Θ.zeta_zw/Θ.Rhostar)) - Θ.gammstar/Θ.bet*(1+Θ.Rhostar)
    Θ.zeta_nn = Θ.gammstar/Θ.bet + Θ.gammstar*Θ.Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Θ.Rhostar)*Θ.muestar*Θ.Gstar*Θ.zeta_Gw/Θ.zeta_zw/Θ.Rhostar
    Θ.zeta_nmue = Θ.gammstar*Θ.Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Θ.Rhostar)*Θ.muestar*Θ.Gstar*(1 - Θ.zeta_Gw*Θ.zeta_zmue/Θ.zeta_zw)
    Θ.zeta_nsigw = Θ.gammstar*Θ.Rkstar/Θ.pistar/exp(Θ.zstar)*(1+Θ.Rhostar)*Θ.muestar*Θ.Gstar*(Θ.zeta_Gsigw-Θ.zeta_Gw/Θ.zeta_zw*Θ.zeta_zsigw)
    
    return Θ
end
