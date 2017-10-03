"""
`init_subspec!(m::SmetsWouters)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::SmetsWouters)

    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss2"
        return ss2!(m)
    else
        error("This subspec is not defined.")
    end
end

"""
```
ss2!(m::SmetsWouters)
```

Initializes subspec 2 of `SmetsWouters`. The bounds for beta-distributed parameters have been changed from
(1e-5, 0.99) to (0.0, 1.0).
"""
function ss2!(m::SmetsWouters)
    m <= parameter(:ζ_p, 0.7813, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1), fixed = false,
                   description = "ζ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label = "\\zeta_p")

    m <= parameter(:ι_p, 0.1865, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.15), fixed = false,
                   description = "ι_p: The weight attributed to last period's inflation in price indexation. (1-ι_p) is the weight attributed to steady-state inflation.",
                   tex_label = "\\iota_p")

    m <= parameter(:h, 0.7205, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.7, 0.1), fixed = false,
                   description = "h: Consumption habit persistence.",
                   tex_label = "h")

    m <= parameter(:ppsi, 0.2648, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.15), fixed = false,
                   description = "ppsi: Utilization costs.",
                   tex_label = "\\psi")

    m <= parameter(:ζ_w, 0.7937, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1), fixed = false,
                   description = "ζ_w: (1-ζ_w) is the probability with which households can freely choose wages in each period. With probability ζ_w, wages increase at a geometrically weighted average of the steady state rate of wage increases and last period's productivity times last period's inflation.",
                   tex_label = "\\zeta_w")

    m <= parameter(:ι_w, 0.4425, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.15), fixed = false,
                   description = "ι_w: The weight attributed to last period's wage in wage indexation. (1-ι_w) is the weight attributed to steady-state wages.",
                   tex_label = "\\iota_w")

    m <= parameter(:ρ, 0.8258, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.75, 0.10), fixed = false,
                   description = "ρ: The degree of inertia in the monetary policy rule.",
                   tex_label = "\\rho_R")

    # Financial frictions parameters

    m <= parameter(:ρ_g, 0.9930, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g")

    m <= parameter(:ρ_b, 0.2703, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   tex_label = "\\rho_b")

    m <= parameter(:ρ_μ, 0.5724, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label = "\\rho_{\\mu}")

    m <= parameter(:ρ_z, 0.9676, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_z: AR(1) coefficient in the technology process.",
                   tex_label = "\\rho_z")

    m <= parameter(:ρ_λ_f, 0.8692, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_λ_f: AR(1) coefficient in the price mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_f}")

    m <= parameter(:ρ_λ_w, 0.9546, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_w}")

    m <= parameter(:ρ_rm, 0.3000, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description = "ρ_rm: AR(1) coefficient in the monetary policy shock process.",
                   tex_label = "\\rho_{r^m}")

    m <= parameter(:η_gz, 0.0500, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.50, 0.20), fixed = false,
                   description = "η_gz: Correlate g and z shocks.",
                   tex_label = "\\eta_{gz}")

    m <= parameter(:η_λ_f, 0.7652, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.50, 0.20), fixed = false,
                   tex_label = "\\eta_{\\lambda_f}")

    m <= parameter(:η_λ_w, 0.8936, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.50, 0.20), fixed = false,
                   description = "η_λ_w: AR(2) coefficient on wage markup shock process.",
                   tex_label = "\\eta_{\\lambda_w}")
end
