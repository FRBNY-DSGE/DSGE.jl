"""
`init_subspec!(m::HetDSGEGovDebt)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::HetDSGEGovDebt)
    if subspec(m) == "ss0"
        return
    # All but g_shock
    elseif subspec(m) == "ss1"
        return ss1!(m)
    # All but b_shock
    elseif subspec(m) == "ss2"
        return ss2!(m)
    # All but μ shock
    elseif subspec(m) == "ss3"
        return ss3!(m)
    # All but z shock
    elseif subspec(m) == "ss4"
        return ss4!(m)
    # All but λ_f shock
    elseif subspec(m) == "ss5"
        return ss5!(m)
    # All but λ_w shock
    elseif subspec(m) == "ss6"
        return ss6!(m)
    # All but rm_shock
    elseif subspec(m) == "ss7"
        return ss7!(m)
    # Fix all steady state relevant parameters
    # in order to estimate non-steady state parameters
    # γ = 0.0
    elseif subspec(m) == "ss8"
        return ss8!(m)
    # Subspec 0, except r has higher prior mean
    elseif subspec(m) == "ss9"
        return ss9!(m)
    # Subspec 8, but with γ = 0.5
    elseif subspec(m) == "ss10"
        return ss10!(m)
    # Do not calibrate for s_H / s_L and zlo
    elseif subspec(m) == "ss11"
        return ss11!(m)
    # Only estimate shocks
    elseif subspec(m) == "ss12"
        return ss12!(m)
    else
        error("This subspec should be a 0")
    end
end

function fix_all_except_sigmas!(m::HetDSGEGovDebt)
    m <= parameter(:α, 0.3, fixed = true, (1e-5, 0.999), (1e-5, 0.999),
                   SquareRoot(), Normal(0.30, 0.05),
                   description = "α: Capital elasticity in the intermediate goods" *
                   "sector's production function (also known as the capital share).",
                   tex_label = "\\alpha")
    # Check this: Previously the description was "Aggregate hours worked"
    m <= parameter(:H, 1.0, (-1000., 1000.), (-1e3, 1e3), Untransformed(),
                   Normal(-45., 5.), fixed = true,
                   description = "Lmean: Mean level of hours", tex_label = "H")
    m <= parameter(:δ, 0.03, fixed = true,
                   description = "δ: The capital depreciation rate", tex_label = "\\delta")
    m <= parameter(:μ_sp, 0.0, fixed = true,
                   description = "μ_sp: The trend in the skill process",
                   tex_label = "\\mu_{sp}")

    m <= parameter(:γ, 0.5, (-5.0, 5.0), (-5., 5.), Untransformed(),
                       Normal(0.4, 0.1), fixed = true, scaling = x -> x/100,
                       description = "γ: The log of the steady-state growth rate of technology",
                       tex_label="100\\gamma")

    m <= parameter(:r, 0.5, (1e-5, 10.0), (1e-5, 10.0), ModelConstructors.Exponential(),
                   GammaAlt(0.5, 0.5), fixed = true, scaling = x -> x/100,
                   description= "r: Quarterly steady-state real interest rate.",
                   tex_label= "100*r^{HetDSGE}")

    m <= parameter(:g, 1/(1-0.01), fixed = true,
                   description = "g_star: 1 - (c_star + i_star)/y_star",
                   tex_label = "g_*")

    m <= parameter(:β_save, 0.0, fixed = true,
                   description = "saving the betas per particle",
                   tex_label = "\\beta_save")
    m <= parameter(:sH_over_sL, 6.33333, (3.0, 9.0), (3.0, 9.0), Untransformed(),
                   Uniform(3.0, 9.0), fixed = true,
                   description = "Ratio of high to low earners", tex_label = "s_H / s_L")

    m <= parameter(:pLH, 0.005, (0.0025, 0.095), (0.0025, 0.095), Untransformed(),
                   Uniform(0.005, 0.095), fixed = true,
                   description = "Prob of going from low to high persistent skill",
                   tex_label = "p(s_L \\mid s_H)")

    m <= parameter(:pHL, 0.03, (0.0025, 0.095), (0.0025, 0.095), Untransformed(),
                   Uniform(0.0025, 0.095), fixed = true,
                   description = "Prob of going from high to low persistent skill",
                   tex_label = "p(s_H \\mid s_L)")

    m <= parameter(:BoverY, 0.26, fixed = true, description = "B / Y", tex_label = "B / Y")

    m <= parameter(:zlo, 0.0323232, (1e-18, 0.8-eps()), (1e-18, 0.8-eps()), Untransformed(),
                   Uniform(1e-18, 0.8-eps()), fixed = true,
                   description = "Lower bound on second income shock to mollify actual income",
                   tex_label = "\\underbar{z}")

    m <= parameter(:zhi, 2-m[:zlo].value, fixed = true,
                   description = "Upper bound on second income shock to mollify actual income",
                   tex_label = "\\bar{z}")

    m <= parameter(:mpc, 0.23395,  fixed = true, tex_label = "MPC")
    m <= parameter(:pc0, 0.071893, fixed = true, description = "Number of people at 0 income",
                   tex_label = "pc0")

    # Give model new parameters
    m <= parameter(:varlinc, 0.0, fixed = true, tex_label = "varlinc",
                   description = "var(log(annual income))")
    m <= parameter(:vardlinc, 0.0, fixed = true, tex_label = "vardlinc",
                   description = "var(log(deviations in annual income))")

    ####################################################
    # Parameters that affect dynamics (not steady-state)
    ####################################################
    m <= parameter(:spp, 4., (-15., 15.), (-15., 15.), Untransformed(),
                   Normal(4., 1.5), fixed = true,
                   description = "S'': The second derivative of households' cost " *
                   "of adjusting investment.", tex_label = "S''")

    m <= parameter(:ϕh, 2., (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(),
                   Normal(2, 0.75), fixed = true, description = "inverse frisch elasticity",
                   tex_label = "\\phi_h")

    m <= parameter(:κ_p, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3),
                   fixed = true, description = "κ_p : The slope of the Price Phillips curve",
                   tex_label = "\\kappa_p")
    m <= parameter(:κ_w, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3),
                   fixed = true, description = "κ_w: The slope of the Wage Phillips curve",
                   tex_label = "\\kappa_w")

    m <= parameter(:ρR , 0.75, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.75, 0.10), fixed = true,
                   description = "ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ψπ , 1.5, (1e-5, 10.), (1e-5, 10.0), ModelConstructors.Exponential(),
                   Normal(1.5, 0.25), fixed = true,
                   description = "ψ1: Weight on inflation gap in monetary policy rule.",
                   tex_label = "\\psi_1")
    m <= parameter(:ψy , 0.5, (-0.5, 0.5), (-0.5, 0.5), Untransformed(),
                   Normal(0.12, 0.05), fixed = true,
                   description = "ψy: Weight on output growth in monetary policy rule",
                   tex_label = "\\psi_y")

    m <= parameter(:δb, 1., (0.0, 1.0), (0.0, 1.0), Untransformed(),
                   Uniform(0.0, 1.0), fixed = true, description = "=1 means balanced budget",
                   tex_label = "\\delta_b")

    # Exogenous processes - autocorrelation
    m <= parameter(:ρ_g, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description = "ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g")
    m <= parameter(:ρ_b, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description = "ρ_b: AR(1) coefficient in intertemporal preference " *
                   "shift process.", tex_label = "\\rho_b")
    m <= parameter(:ρ_μ, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description = "ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label = "\\rho_{\\mu}")
    m <= parameter(:ρ_z, 0.5,(1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description = "ρ_z: AR(1) coefficient in the technology process.",
                   tex_label = "\\rho_z")
    m <= parameter(:ρ_λ_f, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description = "ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_f}")
    m <= parameter(:ρ_λ_w, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description = "ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_w}")
    m <= parameter(:ρ_rm, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = true,
                   description = "ρ_rm: AR(1) coefficient in the monetary policy shock process.",
                   tex_label = "\\rho_{r^m}")

    m <= parameter(:π_star, 0.7000, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(),
                   GammaAlt(0.62, 0.1), fixed = true, scaling = x -> 1 + x/100,
                   description="π_star: steady-state rate of inflation.",
                   tex_label="\\pi_*")

    m <= parameter(:Lmean, -45.9364, (-1000., 1000.), (-1e3, 1e3), Untransformed(),
                   Normal(-45., 5.), fixed = true,
                   description="Lmean: Mean level of hours.", tex_label="\\bar{L}")


end

"""
```
fix_all_but_shocks!(m::HetDSGEGovDebt)
```
Fixes all parameters except the ρ's and the σ's.
"""
function fix_all_but_shocks!(m::HetDSGEGovDebt)
    m <= parameter(:α, 0.3, fixed = true, (1e-5, 0.999), (1e-5, 0.999),
                   SquareRoot(), Normal(0.30, 0.05),
                   description = "α: Capital elasticity in the intermediate goods" *
                   "sector's production function (also known as the capital share).",
                   tex_label = "\\alpha")
    # Check this: Previously the description was "Aggregate hours worked"
    m <= parameter(:H, 1.0, (-1000., 1000.), (-1e3, 1e3), Untransformed(),
                   Normal(-45., 5.), fixed = true,
                   description = "Lmean: Mean level of hours", tex_label = "H")
    m <= parameter(:δ, 0.03, fixed = true,
                   description = "δ: The capital depreciation rate", tex_label = "\\delta")
    m <= parameter(:μ_sp, 0.0, fixed = true,
                   description = "μ_sp: The trend in the skill process",
                   tex_label = "\\mu_{sp}")

    # Exogenous processes - level
    # Uncomment scaling once adjusted properly in the code
    m <= parameter(:γ, 0.5, (-5.0, 5.0), (-5., 5.), Untransformed(),
                       Normal(0.4, 0.1), fixed = true, scaling = x -> x/100,
                       description = "γ: The log of the steady-state growth rate of technology",
                       tex_label="100\\gamma")

    m <= parameter(:r, 0.5, (1e-5, 10.0), (1e-5, 10.0), ModelConstructors.Exponential(),
                   GammaAlt(0.5, 0.5), fixed = true, scaling = x -> x/100,
                   description= "r: Quarterly steady-state real interest rate.",
                   tex_label= "100*r^{HetDSGE}")

    m <= parameter(:g, 1/(1-0.01), fixed = true,
                   description = "g_star: 1 - (c_star + i_star)/y_star",
                   tex_label = "g_*")

    m <= parameter(:β_save, 0.0, fixed = true,
                   description = "saving the betas per particle",
                   tex_label = "\\beta_save")
    m <= parameter(:sH_over_sL, 6.33333, (3.0, 9.0), (3.0, 9.0), Untransformed(),
                   Uniform(3.0, 9.0), fixed = true,
                   description = "Ratio of high to low earners", tex_label = "s_H / s_L")

    m <= parameter(:pLH, 0.005, (0.0025, 0.095), (0.0025, 0.095), Untransformed(),
                   Uniform(0.005, 0.095), fixed = true,
                   description = "Prob of going from low to high persistent skill",
                   tex_label = "p(s_L \\mid s_H)")
    m <= parameter(:pHL, 0.03, (0.0025, 0.095), (0.0025, 0.095), Untransformed(),
                   Uniform(0.0025, 0.095), fixed = true,
                   description = "Prob of going from high to low persistent skill",
                   tex_label = "p(s_H \\mid s_L)")

    m <= parameter(:BoverY, 0.26, fixed = true, description = "B / Y", tex_label = "B / Y")

    m <= parameter(:zlo, 0.0323232, (1e-18, 0.8-eps()), (1e-18, 0.8-eps()), Untransformed(),
                   Uniform(1e-18, 0.8-eps()), fixed = true,
                   description = "Lower bound on second income shock to mollify actual income",
                   tex_label = "\\underbar{z}")

    m <= parameter(:zhi, 2-m[:zlo].value, fixed = true,
                   description = "Upper bound on second income shock to mollify actual income",
                   tex_label = "\\bar{z}")

    m <= parameter(:mpc, 0.23395,  fixed = true, tex_label = "MPC")
    m <= parameter(:pc0, 0.071893, fixed = true, description = "Number of people at 0 income",
                   tex_label = "pc0")

    # Give model new parameters
    m <= parameter(:varlinc, 0.0, fixed = true, tex_label = "varlinc",
                   description = "var(log(annual income))")
    m <= parameter(:vardlinc, 0.0, fixed = true, tex_label = "vardlinc",
                   description = "var(log(deviations in annual income))")

    # Not in m1002
    m <= parameter(:η, 0.0, description = "η: Borrowing constraint (normalized by TFP)",
                   tex_label = "\\eta")

    ####################################################
    # Parameters that affect dynamics (not steady-state)
    ####################################################
    m <= parameter(:spp, 4., (-15., 15.), (-15., 15.), Untransformed(),
                   Normal(4., 1.5), fixed = true,
                   description = "S'': The second derivative of households' cost " *
                   "of adjusting investment.", tex_label = "S''")

    m <= parameter(:ϕh, 2., (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(),
                   Normal(2, 0.75), fixed = true, description = "inverse frisch elasticity",
                   tex_label = "\\phi_h")

    m <= parameter(:κ_p, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3),
                   fixed = true, description = "κ_p : The slope of the Price Phillips curve",
                   tex_label = "\\kappa_p")
    m <= parameter(:κ_w, 0.5, (1e-5, 5.), (1e-5, 5.), SquareRoot(), GammaAlt(0.5, 0.3),
                   fixed = true, description = "κ_w: The slope of the Wage Phillips curve",
                   tex_label = "\\kappa_w")

    m <= parameter(:ρR , 0.75, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
                   BetaAlt(0.75, 0.10), fixed = true,
                   description = "ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")
    m <= parameter(:ψπ , 1.5, (1e-5, 10.), (1e-5, 10.0), ModelConstructors.Exponential(),
                   Normal(1.5, 0.25), fixed = true,
                   description = "ψ1: Weight on inflation gap in monetary policy rule.",
                   tex_label = "\\psi_1")
    m <= parameter(:ψy , 0.5, (-0.5, 0.5), (-0.5, 0.5), Untransformed(),
                   Normal(0.12, 0.05), fixed = true,
                   description = "ψy: Weight on output growth in monetary policy rule",
                   tex_label = "\\psi_y")

    m <= parameter(:δb, 1., (0.0, 1.0), (0.0, 1.0), Untransformed(),
                   Uniform(0.0, 1.0), fixed = true, description = "=1 means balanced budget",
                   tex_label = "\\delta_b")

    # Exogenous processes - autocorrelation
    m <= parameter(:ρ_g, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g")
    m <= parameter(:ρ_b, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_b: AR(1) coefficient in intertemporal preference " *
                   "shift process.", tex_label = "\\rho_B")
    m <= parameter(:ρ_μ, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label = "\\rho_{\\mu}")
    m <= parameter(:ρ_z, 0.5,(1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_z: AR(1) coefficient in the technology process.",
                   tex_label = "\\rho_z")
    m <= parameter(:ρ_λ_f, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_f}")
    m <= parameter(:ρ_λ_w, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_w}")
    m <= parameter(:ρ_rm, 0.5, (1e-5, 1 - 1e-5), (1e-5, 1-1e-5), SquareRoot(),
                   BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_rm: AR(1) coefficient in the monetary policy shock process.",
                   tex_label = "\\rho_{r^m}")

    m <= parameter(:σ_g, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_g: standard dev. of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_b: standard dev. of the intertemporal preference " *
                   "shifter process.", tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_μ: standard dev. of the exogenous marginal efficiency" *
                   " of investment shock process.", tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.15,(1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_z: standard dev. of the process describing the " *
                   "stationary component of productivity.", tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_λ_f: mean of the process that generates the price " *
                   "elasticity of the composite good. Specifically, the elasticity is " *
                   "(1+λ_{f,t})/(λ_{f,t}).", tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_λ_w", tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = false,
                   description = "σ_r_m: standard dev. of the monetary policy shock.",
                   tex_label = "\\sigma_{r^m}")

    m <= parameter(:π_star, 0.7000, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(),
                   GammaAlt(0.62, 0.1), fixed=true, scaling = x -> 1 + x/100,
                   description="π_star: steady-state rate of inflation.",
                   tex_label="\\pi_*")


    m <= parameter(:Lmean, -45.9364, (-1000., 1000.), (-1e3, 1e3), Untransformed(),
                   Normal(-45., 5.), fixed=true,
                   description="Lmean: Mean level of hours.", tex_label="\\bar{L}")

    m <= parameter(:e_y, 0.0, fixed = true,
                   description = "e_y: Measurement error on GDP", tex_label = "e_y")
    m <= parameter(:e_L, 0.0, fixed = true,
                   description = "e_L: Measurement error on hours worked", tex_label = "e_L")
    m <= parameter(:e_w, 0.0, fixed = true, description = "e_w: Measurement error on wages",
                   tex_label ="e_w")
    m <= parameter(:e_π, 0.0, fixed = true,
                   description = "e_π: Measurement error on GDP deflator", tex_label = "e_π")
    m <= parameter(:e_R, 0.0, fixed = true,
                   description = "e_R: Measurement error on nominal rate of interest",
                   tex_label = "e_R")
    m <= parameter(:e_c, 0.0, fixed = true,
                   description = "e_c: Measurement error on consumption", tex_label = "e_c")
    m <= parameter(:e_i, 0.0, fixed = true,
                   description = "e_i: Measurement error on investment", tex_label = "e_i")
end

"""
```
ss1!(m::HetDSGEGovDebt)
```

Initializes subspec 1 of `HetDSGEGovDebt`.
This shuts down all shocks except the government spending process (g shock).
"""
function ss1!(m::HetDSGEGovDebt)
    fix_all_except_sigmas!(m)

    m <= parameter(:σ_b, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.0,(1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f,t}).",
                   tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_w",
                   tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_r_m: The standard deviation of the monetary policy shock.", tex_label = "\\sigma_{r^m}")

end

"""
```
ss2!(m::HetDSGEGovDebt)
```

Initializes subspec 2 of `HetDSGEGovDebt`.
This shuts down all shocks except the preference shifter process (b shock).
"""
function ss2!(m::HetDSGEGovDebt)
    fix_all_except_sigmas!(m)

    m <= parameter(:σ_g, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_μ, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.0,(1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f,t}).",
                   tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_w",
                   tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_r_m: The standard deviation of the monetary policy shock.", tex_label = "\\sigma_{r^m}")

end

"""
```
ss3!(m::HetDSGEGovDebt)
```

Initializes subspec 3 of `HetDSGEGovDebt`.
This shuts down all shocks except the capital adjustment cost process (μ shock).
"""
function ss3!(m::HetDSGEGovDebt)
    fix_all_except_sigmas!(m)

    m <= parameter(:σ_g, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    m <= parameter(:σ_z, 0.0,(1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f,t}).",
                   tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_w",
                   tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_r_m: The standard deviation of the monetary policy shock.", tex_label = "\\sigma_{r^m}")

end

"""
```
ss4(m::HetDSGEGovDebt)
```

Initializes subspec 4 of `HetDSGEGovDebt`.
This shuts down all shocks except the technology process (z shock).
"""
function ss4!(m::HetDSGEGovDebt)

    fix_all_except_sigmas!(m)

    m <= parameter(:σ_g, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_λ_f, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f,t}).",
                   tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_w",
                   tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_r_m: The standard deviation of the monetary policy shock.", tex_label = "\\sigma_{r^m}")

end

"""
```
ss5!(m::HetDSGEGovDebt)
```

Initializes subspec 5 of `HetDSGEGovDebt`.
This shuts down all shocks except the price mark-up shock process (λ_f shock).
"""
function ss5!(m::HetDSGEGovDebt)

    fix_all_except_sigmas!(m)

    m <= parameter(:σ_g, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.0,(1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_w, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_w",
                   tex_label = "\\sigma_{\\lambda_w}")
    m <= parameter(:σ_rm, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_r_m: The standard deviation of the monetary policy shock.", tex_label = "\\sigma_{r^m}")

end

"""
```
ss6!(m::HetDSGEGovDebt)
```

Initializes subspec 6 of `HetDSGEGovDebt`.
This shuts down all shocks except the wage mark-up shock process (λ_w shock).
"""
function ss6!(m::HetDSGEGovDebt)

    fix_all_except_sigmas!(m)

    m <= parameter(:σ_g, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.0,(1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f,t}).",
                   tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_rm, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_r_m: The standard deviation of the monetary policy shock.", tex_label = "\\sigma_{r^m}")

end

"""
```
ss7!(m::HetDSGEGovDebt)
```

Initializes subspec 7 of `HetDSGEGovDebt`.
This shuts down all shocks except the monetary policy shock process (rm shock).
"""
function ss7!(m::HetDSGEGovDebt)

    fix_all_except_sigmas!(m)

    m <= parameter(:σ_g, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")
    m <= parameter(:σ_b, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    m <= parameter(:σ_μ, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")
    m <= parameter(:σ_z, 0.0,(1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label = "\\sigma_z")
    m <= parameter(:σ_λ_f, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f,t}).",
                   tex_label = "\\sigma_{\\lambda_f}")
    m <= parameter(:σ_λ_w, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(),
                   RootInverseGamma(2, 0.10), fixed = true,
                   description = "σ_λ_w",
                   tex_label = "\\sigma_{\\lambda_w}")

end

"""
```
ss8!(m::HetDSGEGovDebt)
```

Initializes subspec 8 of `HetDSGEGovDebt`.
Fixes all parameters relevant to calculating the steady state,
so we can estimate HetDSGEGovDebt at a fixed steady state.
"""
function ss8!(m::HetDSGEGovDebt)
    # Even though some of these are already fixed (and thus redundant
    # to include in this subspec) all are included for the sake of
    # being as explicit as possible about which parameters affect
    # the steady state and thus should be fixed.
    m <= parameter(:α, 0.3, fixed = true,
                   description = "α: Capital elasticity in the intermediate goods" *
                   "sector's production function (also known as the capital share).",
                   tex_label = "\\alpha")
    m <= parameter(:H, 1.0, fixed = true,
                   description = "Aggregate hours worked", tex_label = "H")
    m <= parameter(:δ, 0.03, fixed = true,
                   description = "δ: The capital depreciation rate", tex_label = "\\delta")
    m <= parameter(:μ_sp, 0.0, fixed = true,
                   description = "μ_sp: The trend in the skill process",
                   tex_label = "\\mu_{sp}")

    m <= parameter(:γ, 0.0, fixed = true, scaling = x -> x/100,
                   description = "γ: The log of the steady-state growth rate of technology",
                   tex_label="100\\gamma")

    m <= parameter(:r, 0.5, fixed = true, scaling = x -> x/100,
                   description= "r: Quarterly steady-state real interest rate.",
                   tex_label= "100*r^{HetDSGE}")
    m <= parameter(:g, 1/(1-0.01), fixed = true,
                   description = "g_star: 1 - (c_star + i_star)/y_star",
                   tex_label = "g_*")
    m <= parameter(:sH_over_sL, 6.33333, fixed = true,
                   description = "Ratio of high to low earners", tex_label = "s_H / s_L")
    m <= parameter(:pLH, 0.005, fixed = true,
                   description = "Prob of going from low to high persistent skill",
                   tex_label = "p(s_L \\mid s_H)")
    m <= parameter(:pHL, 0.03, fixed = true,
                   description = "Prob of going from high to low persistent skill",
                   tex_label = "p(s_H \\mid s_L)")

    m <= parameter(:BoverY, 0.26, fixed = true, description = "???", tex_label = "B / Y")

    # Not in m1002
    m <= parameter(:η, 0.0, fixed = true,
                   description = "η: Borrowing constraint (normalized by TFP)",
                   tex_label = "\\eta")

    m <= Setting(:estimate_only_non_steady_state_parameters, true)
end

"""
```
ss9!(m::HetDSGEGovDebt)
```

Initializes subspec 9 of `HetDSGEGovDebt`.
Is subspec 0 except with higher prior mean on r.
"""
function ss9!(m::HetDSGEGovDebt)
    m <= parameter(:r, 1.0, (1e-5, 10.0), (1e-5, 10.0), ModelConstructors.Exponential(),
                   GammaAlt(1.0, 0.5), fixed = false, scaling = x -> x/100,
                   description= "r: Quarterly steady-state real interest rate.",
                   tex_label= "100*r^{HetDSGE}")
end


"""
```
ss10!(m::HetDSGEGovDebt)
```

Initializes subspec 10 of `HetDSGEGovDebt`.
Is subspec 8 except with γ = 0.5.
"""
function ss10!(m::HetDSGEGovDebt)
    ss8!(m)
    m <= parameter(:γ, 0.5, (-5.0, 5.0), (-5., 5.), Untransformed(),
                   Normal(0.4, 0.1), fixed = true, scaling = x -> x/100,
                   description = "γ: The log of the steady-state growth rate of technology",
                   tex_label="100\\gamma")
end

"""
```
ss11!(m::HetDSGEGovDebt)
```

Initializes model for run when we don't calibrate for varlinc and vardlinc.
"""
function ss11!(m::HetDSGEGovDebt)
    m <= Setting(:calibrate_income_targets, false, "Calibrate for varlinc and vardlinc")

    # Set targets
    m <= Setting(:targets, [0.16, 0.10, 0.7, 0.23], "Targets for [MPC, pc0, varlinc, vardlinc]")
    m <= Setting(:target_vars, [:mpc, :pc0, :varlinc, :vardlinc],
                 "Symbols of variables we're targeting")
    m <= Setting(:target_σt, [0.2, 0.1, 0.05, 0.05],
                 "Target \\sigma_t for MPC, pc0, varlinc, and vardlinc")

    # Give model new parameters
    m <= parameter(:varlinc, 0.0, fixed = true, tex_label = "varlinc",
                   description = "var(log(annual income))")
    m <= parameter(:vardlinc, 0.0, fixed = true, tex_label = "vardlinc",
                   description = "var(log(deviations in annual income))")

    # Since not calibrating, we let zlo and s_H / s_L be free parameters
    m <= parameter(:zlo, 1.035e-8, (1e-18, 0.8-eps()), (1e-18, 0.8-eps()), Untransformed(),
                   Uniform(1e-18, 0.8-eps()), fixed = false,
                   description = "Lower bound on second income shock to mollify actual income",
                   tex_label = "\\underbar{z}")
    m <= parameter(:sH_over_sL, 8.99999, (3.0, 9.0), (3.0, 9.0), Untransformed(),
                   Uniform(3.0, 9.0), fixed = false,
                   description = "Ratio of high to low earners", tex_label = "s_H / s_L")
end

"""
```
ss12!(m::HetDSGEGovDebt)
```

Initializes model for run when we don't calibrate for varlinc and vardlinc, and only estimate shocks.
"""
function ss12!(m::HetDSGEGovDebt)
    fix_all_but_shocks!(m)

    m <= Setting(:calibrate_income_targets, false, "Calibrate for varlinc and vardlinc")

    # Only want time series likelihood
    m <= Setting(:ψ_penalty, 0.0)

    # Set targets
    m <= Setting(:targets, [0.16, 0.10, 0.7, 0.23], "Targets for [MPC, pc0, varlinc, vardlinc]")
    m <= Setting(:target_vars, [:mpc, :pc0, :varlinc, :vardlinc],
                 "Symbols of variables we're targeting")
    m <= Setting(:target_σt, [0.2, 0.1, 0.05, 0.05],
                 "Target \\sigma_t for MPC, pc0, varlinc, and vardlinc")

    # Give model new parameters
    m <= parameter(:varlinc, 0.0, fixed = true, tex_label = "varlinc",
                   description = "var(log(annual income))")
    m <= parameter(:vardlinc, 0.0, fixed = true, tex_label = "vardlinc",
                   description = "var(log(deviations in annual income))")

    # Since not calibrating, we let zlo and s_H / s_L be free parameters
    m <= parameter(:zlo, 1.035e-8, (1e-18, 0.8-eps()), (1e-18, 0.8-eps()), Untransformed(),
                   Uniform(1e-18, 0.8-eps()), fixed = true,
                   description = "Lower bound on second income shock to mollify actual income",
                   tex_label = "\\underbar{z}")
    m <= parameter(:sH_over_sL, 8.99999, (3.0, 9.0), (3.0, 9.0), Untransformed(),
                   Uniform(3.0, 9.0), fixed = true,
                   description = "Ratio of high to low earners", tex_label = "s_H / s_L")
end
