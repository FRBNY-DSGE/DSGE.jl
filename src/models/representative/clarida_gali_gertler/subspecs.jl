function init_subspec!(m::ClaridaGaliGertler)
    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss1"
        return ss1!(m)
    elseif subspec(m) == "ss2"
        return ss2!(m)
    elseif subspec(m)
        error("This subspec should be a 0")
    end
end

# Creates a variant of the Clarida-Gali-Gertler (1999) parametrization
function ss1!(m::ClaridaGaliGertler)
    m <= parameter(:γ_Q, .7731, (-1., 10.), (1e-20, 1.), ModelConstructors.Exponential(), Normal(.5, 0.25), fixed = false,
                   description="γ_Q, steady state quarterly growth rate of technology",
                   tex_label="\\gamma)")

    m <= parameter(:π_star, 1.0975, (-1., 10.), (1e-20, 1.), ModelConstructors.Exponential(), Normal(1., .5), fixed = false,
                   description="π_star: Target inflation rate.",
                   tex_label="\\pi*")

    m <= parameter(:rA, 0.7579, (1e-6, 10.), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(0.5, 0.25), fixed = false,
                   description="rA: Target nominal interest rate.",
                   tex_label="r_*")

    m <= parameter(:κ, 0.1, (-10., 10.), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(0.3, .15), fixed = false,
                   description="κ: Composite parameter in New Keynesian Phillips Curve.",
                   tex_label="\\kappa")

    m <= parameter(:τ, 2.3234, (1e-6, 100.), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(2., 0.5), fixed = false,
                   description="τ: The inverse of the intemporal elasticity of substitution.",
                   tex_label="\\tau")


    m <= parameter(:ψ_1, 1.4266, (1e-6, 10.), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(1.5, 0.25), fixed = false,
                   description="ψ_1: The weight on inflation in the monetary policy rule.",
                   tex_label="\\psi_1")
    m <= parameter(:ψ_2, 0.4954, (1e-6, 10.), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(0.125, 0.1), fixed = false,
                   description="ψ_2: The weight on the output gap in the monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ρ_R, 0.6240, (1e-6, 1. - 1e-6), (1e-20, 1.), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_R: AR(1) coefficient on interest rate.",
                   tex_label="\\rho_R")

    m <= parameter(:ρ_g, 0.9126, (1e-6, 1. - 1e-6), (1e-20, 1.), ModelConstructors.SquareRoot(), BetaAlt(0.8, 0.1), fixed = false,
                   description="ρ_g: AR(1) coefficient on g_t = 1/(1 - ζ_t), where ζ_t is government spending as a fraction of output.",
                   tex_label="\\rho_g")

    m <= parameter(:ρ_z, 0.6103, (1e-6, 1. - 1e-6), (1e-20, 1.), ModelConstructors.SquareRoot(), BetaAlt(0.3, 0.1), fixed = false,
                   description="ρ_z: AR(1) coefficient on shocks to the technology growth rate.",
                   tex_label="\\rho_z")

    m <= parameter(:σ_R, 0.4900, (1e-20, 500.), (1e-20, 1.), ModelConstructors.Exponential(), RootInverseGamma(4, .2), fixed = false,
                   description="σ_R: Standard deviation of shocks to the nominal interest rate.",
                   tex_label="\\sigma_R")

    m <= parameter(:σ_g, 1.4594, (1e-20, 500.), (1e-20, 1.), ModelConstructors.Exponential(), RootInverseGamma(4, .5), fixed = false,
                   description="σ_g: Standard deviation of shocks to the government spending process.",
                   tex_label="\\sigma_g")

    m <= parameter(:σ_z, 0.9247, (1e-20, 500.), (1e-20, 1.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.7), fixed = false,
                   description="σ_z: Standard deviation of shocks to the technology growth rate process.",
                   tex_label="\\sigma_z")

    m <= parameter(:e_y, 0., fixed=true,
                   description="e_y: Measurement error on GDP growth.",
                   tex_label="e_y")

    m <= parameter(:e_π, 0., fixed=true,
                   description="e_π: Measurement error on inflation.",
                   tex_label="e_\\pi")

    m <= parameter(:e_R, 0., fixed=true,
                   description="e_R: Measurement error on the interest rate.",
                   tex_label="e_R")
end

function ss2!(m::ClaridaGaliGertler)
    m <= parameter(:γ_Q, .7731, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), Normal(.5, 0.25), fixed = false,
                   description="γ_Q, steady state quarterly growth rate of technology",
                   tex_label="\\gamma)")

    m <= parameter(:π_star, 1.0975, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), Normal(1., .5), fixed = false,
                   description="π_star: Target inflation rate.",
                   tex_label="\\pi*")

    m <= parameter(:rA, 0.7579, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(0.5, 0.25), fixed = false,
                   description="rA: Target nominal interest rate.",
                   tex_label="r_*")

    m <= parameter(:κ, 0.1, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(0.3, .15), fixed = false,
                   description="κ: Composite parameter in New Keynesian Phillips Curve.",
                   tex_label="\\kappa")

    m <= parameter(:τ, 2.3234, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(2., 0.5), fixed = false,
                   description="τ: The inverse of the intemporal elasticity of substitution.",
                   tex_label="\\tau")


    m <= parameter(:ψ_1, 1.4266, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(1.5, 0.25), fixed = false,
                   description="ψ_1: The weight on inflation in the monetary policy rule.",
                   tex_label="\\psi_1")
    m <= parameter(:ψ_2, 0.4954, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), GammaAlt(0.125, 0.1), fixed = false,
                   description="ψ_2: The weight on the output gap in the monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ρ_R, 0.6240, (1e-20, 1.), (1e-20, 1.), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_R: AR(1) coefficient on interest rate.",
                   tex_label="\\rho_R")

    m <= parameter(:ρ_g, 0.9126, (1e-20, 1.), (1e-20, 1.), ModelConstructors.SquareRoot(), BetaAlt(0.8, 0.1), fixed = false,
                   description="ρ_g: AR(1) coefficient on g_t = 1/(1 - ζ_t), where ζ_t is government spending as a fraction of output.",
                   tex_label="\\rho_g")

    m <= parameter(:ρ_z, 0.6103, (1e-20, 1.), (1e-20, 1.), ModelConstructors.SquareRoot(), BetaAlt(0.3, 0.1), fixed = false,
                   description="ρ_z: AR(1) coefficient on shocks to the technology growth rate.",
                   tex_label="\\rho_z")

    m <= parameter(:σ_R, 0.4900, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), RootInverseGamma(4, .2), fixed = false,
                   description="σ_R: Standard deviation of shocks to the nominal interest rate.",
                   tex_label="\\sigma_R")

    m <= parameter(:σ_g, 1.4594, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), RootInverseGamma(4, .5), fixed = false,
                   description="σ_g: Standard deviation of shocks to the government spending process.",
                   tex_label="\\sigma_g")

    m <= parameter(:σ_z, 0.9247, (1e-20, 1e5), (1e-20, 1.), ModelConstructors.Exponential(), RootInverseGamma(4, 0.7), fixed = false,
                   description="σ_z: Standard deviation of shocks to the technology growth rate process.",
                   tex_label="\\sigma_z")

    m <= parameter(:e_y, 0., fixed=true,
                   description="e_y: Measurement error on GDP growth.",
                   tex_label="e_y")

    m <= parameter(:e_π, 0., fixed=true,
                   description="e_π: Measurement error on inflation.",
                   tex_label="e_\\pi")

    m <= parameter(:e_R, 0., fixed=true,
                   description="e_R: Measurement error on the interest rate.",
                   tex_label="e_R")
end
