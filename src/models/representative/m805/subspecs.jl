"""
`init_subspec!(m::Model805)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::Model805)

    if subspec(m) == "ss0"
        # Use σ_r_m as SD for all anticipated shocks without normalizing
        # See measurement equation
        return
    elseif subspec(m) == "ss1"
        # Normalize s.t. sum of variances of anticipated shocks equals the
        # variance of the contemporaneous shock
        # See measurement equation
        return
    elseif subspec(m) =="ss3"
        # Diffuse prior
        ss3!(m)
    elseif subspec(m) == "ss4"
        # Diffuse prior but keep pi* the same
        ss3!(m)
        ss4!(m)
    else
        error("This subspec is not defined.")
    end

end


# Diffuse prior
function ss3!(m::Model805)

    m <= parameter(:α,      0.1687, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     Normal(0.30, 0.15),         fixed=false,
                   description="α: Capital elasticity in the intermediate goods sector's Cobb-Douglas production function.",
                   tex_label="\\alpha")

    m <= parameter(:ζ_p,   0.7467, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),     Distributions.Uniform(0.0, 1.0),          fixed=false,
                   description="ζ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

    m <= parameter(:ι_p,   0.2684, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(), Distributions.Uniform(0.0, 1.0),         fixed=false,
                   description="ι_p: The weight attributed to last period's inflation in price indexation. (1-ι_p) is the weight attributed to steady-state inflation.",

                   tex_label="\\iota_p")

    m <= parameter(:δ,      0.025,  fixed=true,
                   description="δ: The capital depreciation rate.", tex_label="\\delta" )


    m <= parameter(:Upsilon,  1.000,  (0., 10.),     (1e-5, 0.),      ModelConstructors.Exponential(),    GammaAlt(1., 0.5),          fixed=true,
                   description="Υ: The trend evolution of the price of investment goods relative to consumption goods. Set equal to 1.",
                   tex_label="\\mathcal{\\Upsilon}")

    m <= parameter(:Φ,   1.4933, (1., 10.),     (1.00, 10.00),   ModelConstructors.Exponential(),    Normal(1.25, 0.36),         fixed=false,
                   description="Φ: Fixed costs.",
                   tex_label="\\Phi")

    m <= parameter(:S′′,       3.3543, (-15., 15.),   (-15., 15.),     Untransformed(),  Normal(4., 4.5),            fixed=false,
                   description="S'': The second derivative of households' cost of adjusting investment.", tex_label="S\\prime\\prime")

    m <= parameter(:h,        0.4656, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(), Distributions.Uniform(0.0, 1.0),          fixed=false,
                   description="h: Consumption habit persistence.", tex_label="h")

    m <= parameter(:ppsi,     0.7614, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(),  Distributions.Uniform(0.0, 1.0),         fixed=false,
                   description="ppsi: Utilization costs.", tex_label="ppsi")

    m <= parameter(:ν_l,     1.0647, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    Normal(2, 2.25),            fixed=false,
                   description="ν_l: The coefficient of relative risk aversion on the labor term of households' utility function.", tex_label="\\nu_l")

    m <= parameter(:ζ_w,   0.7922, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(), Distributions.Uniform(0.0, 1.0),          fixed=false,
                   description="ζ_w: (1-ζ_w) is the probability with which households can freely choose wages in each period. With probability ζ_w, wages increase at a geometrically weighted average of the steady state rate of wage increases and last period's productivity times last period's inflation.",
                   tex_label="\\zeta_w")

    m <= parameter(:ι_w,   0.5729, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(), Distributions.Uniform(0.0, 1.0),         fixed=false,
                   description="ι_w: No description available.",
                   tex_label="\\iota_w")

    m <= parameter(:λ_w,      1.5000,                                                                               fixed=true,
                   description="λ_w: The wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label="\\lambda_w")

    m <= parameter(:β, 0.7420, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.25, 0.3),        fixed=false,  scaling = x -> 1/(1 + x/100),
                   description="β: Discount rate.",
                   tex_label="\\beta ")


    m <= parameter(:ψ1,  1.8678, (1e-5, 10.),   (1e-5, 10.00),   ModelConstructors.Exponential(),    Normal(1.5, 0.75),          fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2,  0.0715, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.15),         fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2131, (-0.5, 0.5),   (-0.5, 0.5),     Untransformed(),  Normal(0.12, 0.15),         fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")
 ## The old prior for π* in SW is GammaAlt(0.62, 0.1) with diffuse GammaAlt(0.62, 0.3) and the old prior here was GammaAlt(0.75, 0.4) so we already had a higher variance, but in keeping with the 'multiply by 3' rule, I'm multiplying by 3 anyways.
    m <= parameter(:π_star,   0.6231, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.75, 1.2),        fixed=false,  scaling = x -> 1 + x/100,
                   description="π_star: The steady-state rate of inflation.",
                   tex_label="\\pi_*")

m <= parameter(:σ_c, 1.5073, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    Normal(1.5, 1.11),          fixed=false,
               description="σ_c: No description available.",
               tex_label="\\sigma_{c}")

m <= parameter(:ρ,      .8519, (1e-5, 0.999), (1e-5, 0.999),   SquareRoot(), Distributions.Uniform(0.0, 1.0),        fixed=false,
               description="ρ: The degree of inertia in the monetary policy rule.",
               tex_label="\\rho")

m <= parameter(:ϵ_p,     10.000,                                                                               fixed=true,
               description="ϵ_p: No description available.",
               tex_label="\\varepsilon_{p}")

m <= parameter(:ϵ_w,     10.000,                                                                               fixed=true,
               description="ϵ_w: No description available.",
               tex_label="\\varepsilon_{w}")


# exogenous processes - level
m <= parameter(:γ,      0.3085, (-5.0, 5.0),     (-5., 5.),     Untransformed(), Normal(0.4, 0.3),            fixed=false, scaling = x -> x/100,
               description="γ: The log of the steady-state growth rate of technology.",
               tex_label="\\gamma")

m <= parameter(:Lmean,  -45., (-1000., 1000.), (-1e3, 1e3),   Untransformed(), Normal(-45, 15),   fixed=false,
               description="Lmean: No description available.",
               tex_label="Lmean")

m <= parameter(:g_star,    0.1800,                                                                               fixed=true,
               description="g_star: No description available.",
               tex_label="g_*")

# exogenous processes - autocorrelation
m <= parameter(:ρ_g,      0.9930, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),           fixed=false,
               description="ρ_g: AR(1) coefficient in the government spending process.",
               tex_label="\\rho_g")

m <= parameter(:ρ_b,      0.8100, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),           fixed=false,
               description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
               tex_label="\\rho_b")

m <= parameter(:ρ_μ,     0.7790, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),           fixed=false,
               description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
               tex_label="\\rho_{\\mu}")

m <= parameter(:ρ_z,      0.9676, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),           fixed=false,
               description="ρ_z: AR(1) coefficient in the technology process.",
               tex_label="\\rho_z")

m <= parameter(:ρ_λ_f,    0.8372, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),           fixed=false,
               description="ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
               tex_label="\\rho_{\\lambda_f}")

m <= parameter(:ρ_λ_w,    0.9853, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),           fixed=false,
               description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.", # CHECK THIS
               tex_label="\\rho_{\\lambda_w}")

m <= parameter(:ρ_rm,     0.9000, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),           fixed=false,
               description="ρ_rm: AR(1) coefficient in the monetary policy shock process.", # CHECK THIS
               tex_label="\\rho_{rm}")

m <= parameter(:ρ_π_star,     0.9900, (1e-5, 0.999),   (1e-5, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),
               fixed=true,
               description="ρ_π_star: AR(1) coefficient in the π* process.",
               tex_label="\\rho_{\\pi^*}")
# exogenous processes - standard deviation
m <= parameter(:σ_g,      0.6090, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_g: The standard deviation of the government spending process.",
               tex_label="\\sigma_{g}")

m <= parameter(:σ_b,      0.1818, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_b: No description available.",
               tex_label="\\sigma_{b}")

m <= parameter(:σ_μ,     0.4601, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
               tex_label="\\sigma_{\\mu}")

m <= parameter(:σ_z,      0.4618, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_z: No description available.",
               tex_label="\\sigma_{z}")

m <= parameter(:σ_λ_f,    0.1455, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_λ_f: The mean of the process that generates the price elasticity of the composite good.  Specifically, the elasticity is (1+λ_{f,t})/(λ_{f_t}).",
               tex_label="\\sigma_{\\lambda_f}")

m <= parameter(:σ_λ_w,    0.2089, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_λ_w: No description available.",
               tex_label="\\sigma_{\\lambda_w}")

m <= parameter(:σ_rm,     0.2397, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_rm: No description available.",
               tex_label="\\sigma_{rm}")

m <= parameter(:σ_π_star,      0.0300, (1e-8, 5.),      (1e-8, 5.),    ModelConstructors.Exponential(),   RootInverseGamma(2., 0.10),  fixed=false,
               description="σ_π_star: The standard deviation of the π* process.",
               tex_label="\\sigma_{\\pi^*}")

m <= parameter(:η_gz,       0.5632, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
               Distributions.Uniform(0.0, 1.0), fixed=false,
               description="η_gz: Correlate g and z shocks.",
               tex_label="\\eta_{gz}")

m <= parameter(:η_λ_f,      0.7652, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
               Distributions.Uniform(0.0, 1.0),         fixed=false,
               description="η_λ_f: No description available.",
               tex_label="\\eta_{\\lambda_f}")

m <= parameter(:η_λ_w,      0.8936, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(),
               Distributions.Uniform(0.0, 1.0),         fixed=false,
               description="η_λ_w: AR(2) coefficient on wage markup shock process.",
               tex_label="\\eta_{\\lambda_w}")


end


function ss4!(m::Model805)
    m <= parameter(:π_star,   0.6231, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.75, 0.4),        fixed=false,  scaling = x -> 1 + x/100,
                   description="π_star: The steady-state rate of inflation.",
                   tex_label="\\pi_*")
end
