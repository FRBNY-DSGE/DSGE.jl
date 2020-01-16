"""
`init_subspec!(m::SmetsWoutersOrig)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::SmetsWoutersOrig)

    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss1"
        # Fix γ at 0.4312
        ss1!(m)
    elseif subspec(m) == "ss2"
        # Diffuse prior
        ss2!(m)
    elseif subspec(m) == "ss21"
        return ss21!(m)
    elseif subspec(m) == "ss22"
        return ss22!(m)
    elseif subspec(m) == "ss23"
        return ss23!(m)
    elseif subspec(m) == "ss24"
        return ss24!(m)
    elseif subspec(m) == "ss25"
        return ss25!(m)
    elseif subspec(m) == "ss26"
        return ss26!(m)
    elseif subspec(m) == "ss27"
        return ss27!(m)
    elseif subspec(m) == "ss28"
        return ss28!(m)
    elseif subspec(m) == "ss29"
        return ss29!(m)
        error("This subspec is not defined.")
    end
end

function ss1!(m::SmetsWoutersOrig)
    m <= parameter(:γ, 0.4312, (-5.0, 5.0), (-5., 5.), Untransformed(), Normal(0.4, 0.1),
                   fixed=true, scaling = x -> 1 + x/100,
                   description="γ: The log of the steady-state growth rate of technology.",
                   tex_label="\\gamma")
end

# Diffuse prior
function ss2!(m::SmetsWoutersOrig)
 m <= parameter(:S′′, 6.3325, (-15., 15.), (-15., 15.), Untransformed(), Distributions.Normal(4., 4.5), fixed=false,
                   description="S'': The second derivative of households' cost of adjusting investment.",
                   tex_label="S''")

    m <= parameter(:σ_c, 1.2312, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Distributions.Normal(1.5, 1.11), fixed=false,
                   tex_label="\\sigma_{c}")

    m <= parameter(:h, 0.7205, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
                   description="h: Consumption habit persistence.",
                   tex_label="h")

    m <= parameter(:λ_w, 1.5000, fixed=true,
                   description="λ_w: The wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label="\\lambda_w")

    m <= parameter(:ζ_w, 0.7937, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
                   description="ζ_w: (1-ζ_w) is the probability with which households can freely choose wages in each period. With probability ζ_w, wages increase at a geometrically weighted average of the steady state rate of wage increases and last period's productivity times last period's inflation.",
                   tex_label="\\zeta_w")

    m <= parameter(:ϵ_w, 10.000, fixed=true,
                   description="ϵ_w: Curvature parameter in the Kimball aggregator for wages.",
                   tex_label="\\epsilon_{w}")

    m <= parameter(:ν_l, 2.8401, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Distributions.Normal(2.0, 2.25), fixed=false,
                   description="ν_l: The coefficient of relative risk aversion on the labor term of households' utility function.", tex_label="\\nu_l")

    m <= parameter(:δ, 0.025, fixed=true,
                   description="δ: The capital depreciation rate.",
                   tex_label="\\delta" )

    m <= parameter(:ζ_p, 0.7813, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
                   description="ζ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ϵ_p, 10.000, fixed=true,
                   description="ϵ_p: Curvature parameter in the Kimball aggregator for prices.",
                   tex_label="\\epsilon_{p}")

    m <= parameter(:ι_w, 0.4425, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
                   description="ι_w: The weight attributed to last period's wage in wage indexation. (1-ι_w) is the weight attributed to steady-state wages.",
                   tex_label="\\iota_w")
    m <= parameter(:ι_p, 0.3291, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
                   description="ι_p: The weight attributed to last period's inflation in price indexation. (1-ι_p) is the weight attributed to steady-state inflation.",
                   tex_label="\\iota_p")

    m <= parameter(:ppsi, 0.2648, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
                   description="ppsi: Utilization costs.",
                   tex_label="\\psi")

    m <= parameter(:Φ, 1.4672, (1., 10.), (1.00, 10.00), ModelConstructors.Exponential(), Distributions.Normal(1.25, 0.36), fixed=false,
                   description="Φ: Fixed costs.",
                   tex_label="\\Phi")

    m <= parameter(:ψ1, 1.7985, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Distributions.Normal(1.5, 0.75), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ρ, .8258, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ψ2, 0.0893, (-0.5, 0.5), (-0.5, 0.5), Untransformed(), Distributions.Normal(0.12, 0.15), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2239, (-0.5, 0.5), (-0.5, 0.5), Untransformed(), Distributions.Normal(0.12, 0.15), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:π_star, 0.7000, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), GammaAlt(0.62, 0.3), fixed=false, scaling = x -> 1 + x/100,
                   description="π_star: The steady-state rate of inflation.",
                   tex_label="\\pi_*")

m <= parameter(:β, 0.7420, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), GammaAlt(0.25, 0.3), fixed=false, scaling = x -> 1/(1 + x/100),
               description="β: Discount rate.",
               tex_label="\\beta ")

m <= parameter(:Lmean, -1., (-5., 5.), (-5., 5.), Untransformed(), Distributions.Normal(0.0, 6.0), fixed=false,
               description="Lmean: Mean level of hours.",
               tex_label="\\bar{L}")

m <= parameter(:γ, 0.4312, (-5.0, 5.0), (-5., 5.), Untransformed(), Distributions.Normal(0.4, 0.3), fixed=false, scaling = x -> 1 + x/100,
               description="γ: The log of the steady-state growth rate of technology.",
               tex_label="\\gamma")

m <= parameter(:α, 0.24, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Normal(0.30, 0.15), fixed=false,
               description="α: Capital elasticity in the intermediate goods sector's Cobb-Douglas production function.",
               tex_label="\\alpha")

m <= parameter(:g_star, 0.1800, fixed=true,
               description="g_star: 1 - (c_star + i_star)/y_star.",
               tex_label="g_*")

# exogenous processes - autocorrelation
m <= parameter(:ρ_z, 0.9676, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="ρ_z: AR(1) coefficient in the technology process.",
               tex_label="\\rho_z")

m <= parameter(:ρ_b, 0.2703, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
               tex_label="\\rho_b")
#added 4 9's to upper bound for the following parameter
m <= parameter(:ρ_g, 0.9930, (1e-5, 0.9999999), (1e-5, 0.9999999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="ρ_g: AR(1) coefficient in the government spending process.",
               tex_label="\\rho_g")

m <= parameter(:ρ_μ, 0.5724, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
               tex_label="\\rho_{\\mu}")

m <= parameter(:ρ_rm, 0.3000, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="ρ_rm: AR(1) coefficient in the monetary policy shock process.", # CHECK THIS
               tex_label="\\rho_{r^m}")

m <= parameter(:ρ_λ_f, 0.8692, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
               tex_label="\\rho_{\\lambda_f}")


m <= parameter(:ρ_λ_w, 0.9546, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.", # CHECK THIS
               tex_label="\\rho_{\\lambda_f}")

m <= parameter(:η_λ_f, 0.7652, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="η_λ_f: Moving average component in the price markup shock.",
               tex_label="\\eta_{\\lambda_f}")

m <= parameter(:η_λ_w, 0.8936, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="η_λ_w: Moving average component in the wage markup shock.",
               tex_label="\\eta_{\\lambda_w}")

m <= parameter(:η_gz, 0.0500, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0), fixed=false,
               description="η_gz: Correlate g and z shocks.",
               tex_label="\\eta_{gz}")

# exogenous processes - standard deviation
m <= parameter(:σ_z, 0.4618, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
               description="σ_z: The standard deviation of the process describing the stationary component of productivity.",
               tex_label="\\sigma_{z}")

m <= parameter(:σ_b, 0.1818, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
               description="σ_b: The standard deviation of the intertemporal preference shifter process.",
               tex_label="\\sigma_{b}")

m <= parameter(:σ_g, 0.6090, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
               description="σ_g: The standard deviation of the government spending process.",
               tex_label="\\sigma_{g}")

m <= parameter(:σ_μ, 0.4601, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
               description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
               tex_label="\\sigma_{\\mu}")

m <= parameter(:σ_rm, 0.2397, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
               description="σ_r_m: The standard deviation of the monetary policy shock.",
               tex_label="\\sigma_{r^m}")


m <= parameter(:σ_λ_f, 0.1455, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
               description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
               tex_label="\\sigma_{\\lambda_f}")

m <= parameter(:σ_λ_w, 0.2089, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
               tex_label="\\sigma_{\\lambda_w}")

end

"""
```
ss21!(m::SmetsWoutersOrig)
```

Has two ζ_p parameters for the regime-switching estimation of DSGE. Loose prior!
"""
function ss21!(m::SmetsWoutersOrig)
    m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.3), fixed=false,
                   description="ζ_p: The Calvo parameter (Regime 1). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ζ_p2, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.3), fixed=false,
                   description="ζ_p2: The Calvo parameter (Regime 2). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p^2")
end

"""
```
ss22!(m::SmetsWoutersOrig)
```

Has two ζ_p parameters for the regime-switching estimation of DSGE. Standard prior!
"""
function ss22!(m::SmetsWoutersOrig)
    m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p: The Calvo parameter (Regime 1). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ζ_p2, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p2: The Calvo parameter (Regime 2). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p^2")
end

"""
```
ss23!(m::SmetsWoutersOrig)
```

Has two of the MP rule parameters for the regime-switching estimation of DSGE. Loose prior!
"""
function ss23!(m::SmetsWoutersOrig)
    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.75), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.30), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ψ12, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.75), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ22, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ32, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ2, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.30), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")
end


"""
```
ss24!(m::SmetsWoutersOrig)
```

Has two of the MP rule parameters for the regime-switching estimation of DSGE. Standard prior
"""
function ss24!(m::SmetsWoutersOrig)
    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ψ12, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ22, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ32, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ2, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")
end

"""
```
ss25!(m::SmetsWoutersOrig)
```
Has two of the MP rule parameters and zeta_p for the regime-switching estimation of DSGE. Loose prior!
"""
function ss25!(m::SmetsWoutersOrig)
    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.75), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.30), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ψ12, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.75), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ22, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ32, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ2, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.30), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ζ_p2, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.3), fixed=false,
                   description="ζ_p: The Calvo parameter (Regime 1). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")
    m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.3), fixed=false,
                   description="ζ_p2: The Calvo parameter (Regime 2). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p^2")
end

"""
```
ss26!(m::SmetsWoutersOrig)
```

Has two of the MP rule parameters and zeta_p for the regime-switching estimation of DSGE. Standard prior
"""
function ss26!(m::SmetsWoutersOrig)
    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ψ12, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ22, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ32, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ2, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p: The Calvo parameter (Regime 1). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")


    m <= parameter(:ζ_p2, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p2: The Calvo parameter (Regime 2). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p^2")
end

"""
```
ss27!(m::SmetsWoutersOrig)
```

Has two of each shock's sigma for the regime-switching esitmation of DSGE. Standard prior (because loose prior is the same for these)
"""
function ss27!(m::SmetsWoutersOrig)
    m <= parameter(:σ_g2, 2.5230, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label="\\sigma_{g}")

    m <= parameter(:σ_b2, 0.0292, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label="\\sigma_{b}")

    m <= parameter(:σ_μ2, 0.4559, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label="\\sigma_{\\mu}")

    m <= parameter(:σ_ztil2, 0.6742, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_ztil: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label="\\sigma_{\\tilde{z}}")

    m <= parameter(:σ_λ_f2, 0.1314, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
                   tex_label="\\sigma_{\\lambda_f}")

    m <= parameter(:σ_λ_w2, 0.3864, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{\\lambda_w}")

    m <= parameter(:σ_r_m2, 0.2380, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_r_m: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{r^m}")

    m <= parameter(:σ_σ_ω2, 0.0428, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=false,
                   description="σ_σ_ω: The standard deviation of entrepreneurs' capital productivity follows an exogenous process with standard deviation σ_σ_ω.",
                   tex_label="\\sigma_{\\sigma_\\omega}")

    m <= parameter(:σ_μ_e2, 0.0000, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=true,
                   description="σ_μ_e: Exogenous bankrupcy costs follow an exogenous process with standard deviation σ_μ_e.",
                   tex_label="\\sigma_{\\mu_e}")

    m <= parameter(:σ_γ2, 0.0000, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.01), fixed=true,
                   description="σ_γ: The fraction of entrepreneurs surviving period t follows an exogenous process with standard deviation σ_γ.",
                   tex_label="\\sigma_{\\gamma}")

    m <= parameter(:σ_π_star2, 0.0269, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(6, 0.03), fixed=false,
                   description="σ_π_star: The standard deviation of the inflation target.",
                   tex_label="\\sigma_{\\pi_*}")

    m <= parameter(:σ_lr2, 0.1766, (1e-8,10.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.75), fixed=false,
                   tex_label="\\sigma_{10y}")

    m <= parameter(:σ_z_p2, 0.1662, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_z_p: The standard deviation of the shock to the permanent component of productivity.",
                   tex_label="\\sigma_{z^p}")

    m <= parameter(:σ_tfp2, 0.9391, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{tfp}")

    m <= parameter(:σ_gdpdef2, 0.1575, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdpdef}")

    m <= parameter(:σ_corepce2, 0.0999, (1e-8, 5.),(1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{pce}")

    m <= parameter(:σ_gdp2, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdp}")

    m <= parameter(:σ_gdi2, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdi}")
end

"""
```
ss28!(m::SmetsWoutersOrig)
```

 Has two of each shock's sigma and two zeta_p and two of each MP rule's parameters for the regime-switching esitmation of DSGE. Loose prior (but the sigmas have same as normal because their loose prior is same as standard)
"""
function ss28!(m::SmetsWoutersOrig)
    m <= parameter(:σ_g2, 2.5230, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label="\\sigma_{g}")

    m <= parameter(:σ_b2, 0.0292, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label="\\sigma_{b}")

    m <= parameter(:σ_μ2, 0.4559, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label="\\sigma_{\\mu}")

    m <= parameter(:σ_ztil2, 0.6742, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_ztil: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label="\\sigma_{\\tilde{z}}")

    m <= parameter(:σ_λ_f2, 0.1314, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
                   tex_label="\\sigma_{\\lambda_f}")

    m <= parameter(:σ_λ_w2, 0.3864, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{\\lambda_w}")

    m <= parameter(:σ_r_m2, 0.2380, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_r_m: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{r^m}")

    m <= parameter(:σ_σ_ω2, 0.0428, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=false,
                   description="σ_σ_ω: The standard deviation of entrepreneurs' capital productivity follows an exogenous process with standard deviation σ_σ_ω.",
                   tex_label="\\sigma_{\\sigma_\\omega}")

    m <= parameter(:σ_μ_e2, 0.0000, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=true,
                   description="σ_μ_e: Exogenous bankrupcy costs follow an exogenous process with standard deviation σ_μ_e.",
                   tex_label="\\sigma_{\\mu_e}")

    m <= parameter(:σ_γ2, 0.0000, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.01), fixed=true,
                   description="σ_γ: The fraction of entrepreneurs surviving period t follows an exogenous process with standard deviation σ_γ.",
                   tex_label="\\sigma_{\\gamma}")

    m <= parameter(:σ_π_star2, 0.0269, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(6, 0.03), fixed=false,
                   description="σ_π_star: The standard deviation of the inflation target.",
                   tex_label="\\sigma_{\\pi_*}")

    m <= parameter(:σ_lr2, 0.1766, (1e-8,10.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.75), fixed=false,
                   tex_label="\\sigma_{10y}")

    m <= parameter(:σ_z_p2, 0.1662, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_z_p: The standard deviation of the shock to the permanent component of productivity.",
                   tex_label="\\sigma_{z^p}")

    m <= parameter(:σ_tfp2, 0.9391, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{tfp}")

    m <= parameter(:σ_gdpdef2, 0.1575, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdpdef}")

    m <= parameter(:σ_corepce2, 0.0999, (1e-8, 5.),(1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{pce}")

    m <= parameter(:σ_gdp2, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdp}")

    m <= parameter(:σ_gdi2, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdi}")

    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.75), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
               description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
               tex_label="\\psi_3")

m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.30), fixed=false,
               description="ρ: The degree of inertia in the monetary policy rule.",
               tex_label="\\rho_R")

m <= parameter(:ψ12, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.75), fixed=false,
               description="ψ₁: Weight on inflation gap in monetary policy rule.",
               tex_label="\\psi_1")

m <= parameter(:ψ22, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
               description="ψ₂: Weight on output gap in monetary policy rule.",
               tex_label="\\psi_2")

m <= parameter(:ψ32, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.15), fixed=false,
               description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
               tex_label="\\psi_3")

m <= parameter(:ρ2, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.30), fixed=false,
               description="ρ: The degree of inertia in the monetary policy rule.",
               tex_label="\\rho_R")

m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.3), fixed=false,
               description="ζ_p: The Calvo parameter (Regime 1). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
               tex_label="\\zeta_p")

m <= parameter(:ζ_p2, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.3), fixed=false,
               description="ζ_p2: The Calvo parameter (Regime 2). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
               tex_label="\\zeta_p^2")
end

"""
```
ss29!(m::SmetsWoutersOrig)
```

Has two of each shock's sigma and two zeta_p and two of each MP rule's parameters for the regime-switching esitmation of DSGE. Standard prior (but the sigmas have same as normal because their loose prior is same as standard)
"""
function ss29!(m::SmetsWoutersOrig)
    m <= parameter(:σ_g2, 2.5230, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label="\\sigma_{g}")

    m <= parameter(:σ_b2, 0.0292, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label="\\sigma_{b}")

    m <= parameter(:σ_μ2, 0.4559, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label="\\sigma_{\\mu}")

    m <= parameter(:σ_ztil2, 0.6742, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_ztil: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label="\\sigma_{\\tilde{z}}")

    m <= parameter(:σ_λ_f2, 0.1314, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good. Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
                   tex_label="\\sigma_{\\lambda_f}")

    m <= parameter(:σ_λ_w2, 0.3864, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{\\lambda_w}")

    m <= parameter(:σ_r_m2, 0.2380, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_r_m: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{r^m}")

    m <= parameter(:σ_σ_ω2, 0.0428, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=false,
                   description="σ_σ_ω: The standard deviation of entrepreneurs' capital productivity follows an exogenous process with standard deviation σ_σ_ω.",
                   tex_label="\\sigma_{\\sigma_\\omega}")

    m <= parameter(:σ_μ_e2, 0.0000, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.05), fixed=true,
                   description="σ_μ_e: Exogenous bankrupcy costs follow an exogenous process with standard deviation σ_μ_e.",
                   tex_label="\\sigma_{\\mu_e}")

    m <= parameter(:σ_γ2, 0.0000, (1e-7,100.), (1e-5, 0.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.01), fixed=true,
                   description="σ_γ: The fraction of entrepreneurs surviving period t follows an exogenous process with standard deviation σ_γ.",
                   tex_label="\\sigma_{\\gamma}")

    m <= parameter(:σ_π_star2, 0.0269, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(6, 0.03), fixed=false,
                   description="σ_π_star: The standard deviation of the inflation target.",
                   tex_label="\\sigma_{\\pi_*}")

    m <= parameter(:σ_lr2, 0.1766, (1e-8,10.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.75), fixed=false,
                   tex_label="\\sigma_{10y}")

    m <= parameter(:σ_z_p2, 0.1662, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_z_p: The standard deviation of the shock to the permanent component of productivity.",
                   tex_label="\\sigma_{z^p}")

    m <= parameter(:σ_tfp2, 0.9391, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{tfp}")

    m <= parameter(:σ_gdpdef2, 0.1575, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdpdef}")

    m <= parameter(:σ_corepce2, 0.0999, (1e-8, 5.),(1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{pce}")

    m <= parameter(:σ_gdp2, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdp}")

    m <= parameter(:σ_gdi2, 0.1, (1e-8, 5.),(1e-8, 5.),ModelConstructors.Exponential(),RootInverseGamma(2, 0.10), fixed=false,
                   tex_label="\\sigma_{gdi}")

    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")

    m <= parameter(:ψ12, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ22, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ32, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ2, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
               tex_label="\\ρ_R")

    m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p: The Calvo parameter (Regime 1). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ζ_p2, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p2: The Calvo parameter (Regime 2). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p^2")
end
