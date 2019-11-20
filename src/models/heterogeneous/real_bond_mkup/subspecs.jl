"""
`init_subspec!(m::RealBondMkup)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::RealBondMkup)
    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss1"
        return ss1!(m)
    elseif subspec(m) == "ss2"
        return ss2!(m)
    elseif subspec(m) == "ss3"
        return ss3!(m)
    elseif subspec(m) == "ss4"
        return ss4!(m)
    elseif subspec(m) == "ss5"
        return ss5!(m)
    elseif subspec(m) == "ss6"
        return ss6!(m)
    else
        error("This subspec should be a 0")
    end
end

"""
```
ss1!(m::RealBondMkup)
```

Initializes subspec 1 of `RealBondMkup`.
This shuts down all shocks except the technology process (z shock).
"""
function ss1!(m::RealBondMkup)
    m <= parameter(:ρ_z, 0.0, fixed = true,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
    m <= parameter(:ρ_mon, 0.0, fixed = true,
                   description="ρ_mon: Persistence of the monetary policy shock",
                   tex_label="\\rho_{mon}")
    m <= parameter(:ρ_mkp, 0.0, fixed = true,
                   description="ρ_mkp: Persistence of the markup shock",
                   tex_label="\\rho_{mkp}")

    m <= parameter(:σ_mon, 0.0, fixed = true,
                   description="σ_mon: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{mon}")
    m <= parameter(:σ_mkp, 0.0, fixed = true,
                   description="σ_mkp: The standard deviation of the markup shock.",
                   tex_label="\\sigma_{mkp}")
end

"""
```
ss2!(m::RealBondMkup)
```

Initializes subspec 2 of `RealBondMkup`.
This shuts down all shocks except the markup shock process (mkp shock).
"""
function ss2!(m::RealBondMkup)
    m <= parameter(:ρ_z, 0.0, fixed = true,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
    m <= parameter(:ρ_mon, 0.0, fixed = true,
                   description="ρ_mon: Persistence of the monetary policy shock",
                   tex_label="\\rho_{mon}")
    m <= parameter(:ρ_mkp, 0.0, fixed = true,
                   description="ρ_mkp: Persistence of the markup shock",
                   tex_label="\\rho_{mkp}")

    m <= parameter(:σ_z, 0.0, fixed = true,
                   description="σ_z: The standard deviation of the technology process.",
                   tex_label="\\sigma_z")
    m <= parameter(:σ_mon, 0.0, fixed = true,
                   description="σ_mon: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{mon}")
end

"""
```
ss3!(m::RealBondMkup)
```

Initializes subspec 3 of `RealBondMkup`.
This shuts down all shocks except the monetary policy process (mon shock).
"""
function ss3!(m::RealBondMkup)
    m <= parameter(:ρ_z, 0.0, fixed = true,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
    m <= parameter(:ρ_mon, 0.0, fixed = true,
                   description="ρ_mon: Persistence of the monetary policy shock",
                   tex_label="\\rho_z")
    m <= parameter(:ρ_mkp, 0.0, fixed = true,
                   description="ρ_mkp: Persistence of the markup shock",
                   tex_label="\\rho_z")

    m <= parameter(:σ_z, 0.0, fixed = true,
                   description="σ_z: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_z")
    m <= parameter(:σ_mkp, 0.0, fixed = true,
                   description="σ_mkp: The standard deviation of the markup shock.",
                   tex_label="\\sigma_{mkp}")
end


"""
```
ss4!(m::RealBondMkup)
```

Initializes subspec 4 of `RealBondMkup`.
This fixes all parameters except for the sigmas. Used for estimating only the sigmas
"""
function ss4!(m::RealBondMkup)
    # Steady State Parameters
    m <= parameter(:R, 1.04, fixed = true,
                   description = "R: Steady-state gross real interest rate.", tex_label = "R")
    m <= parameter(:γ, 1.0, fixed = true,
                   description = "γ: CRRA Parameter.", tex_label = "\\gamma")
    m <= parameter(:ν, 1.0, fixed = true,
                   description = "Inverse Frisch elasticity of labor supply.", tex_label = "\\nu")
    m <= parameter(:abar, -0.5, fixed = true,
                   description = "Borrowing floor.", tex_label = "\\bar{a}")

    # Ar for observables (non-steadystate)
    m <= parameter(:ρ_z, 0.03081565350294113, (1e-3, 0.999), (1e-3, 0.999), SquareRoot(), BetaAlt(0.5, 0.2),
                   fixed=true,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
    m <= parameter(:ρ_mon, 0.5, (1e-3, 0.999), (1e-3, 0.999), SquareRoot(), BetaAlt(0.5, 0.2),
                   fixed = true,
                   description = "ρ_mon: Persistence of monetary policy shock",
                   tex_label = "\\rho_{mon}")
    m <= parameter(:ρ_mkp, 0.5, (1e-3, 0.999), (1e-3, 0.999), SquareRoot(), BetaAlt(0.5, 0.2),
                   fixed = true,
                   tex_label = "\\rho_{mkp}", description = "ρ_mkp: Persistence of the markup shock")

    # More non-steadystate parameters
    m <= parameter(:ρ_tay, 0.5, (1e-3, 0.999), (1e-3, 0.999), SquareRoot(), BetaAlt(0.5, 0.2),
                   fixed = true,
                   description = "ρ_tay: Persistence in the taylor rule",
                   tex_label = "rho_{tay}")
    m <= parameter(:κ, 1.0, (1e-3, 0.999), (1e-3, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),
                   fixed = true,
                   description = "κ: The slope of the Phillips curve",
                   tex_label = "\\kappa")
    m <= parameter(:phipi, 1.5, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25),
                   fixed = true,
                   description = "phipi: The slope of the taylor rule",
                   tex_label = "\\phi_\\pi")

    # More steadystate parameters
    m <= parameter(:μ_s, 0., fixed = true, description = "μ_s: Mu of log normal in income")
    m <= parameter(:σ_s, 0.1, fixed = true,
                   description = "σ_s: Sigma of log normal in income")
    m <= parameter(:e_y, 1e-3, fixed = true, description = "e_y: Measurement error on GDP",
                   tex_label = "e_y")
end

"""
```
ss5!(m::RealBondMkup)
```

Initializes subspec 5 of `RealBondMkup`.
This fixes all parameters except for the sigmas and rhos. Used for estimating only the sigmas and rhos.
"""
function ss5!(m::RealBondMkup)
    # Steady State Parameters
    m <= parameter(:R, 1.04, fixed = true,
                   description = "R: Steady-state gross real interest rate.", tex_label = "R")
    m <= parameter(:γ, 1.0, fixed = true,
                   description = "γ: CRRA Parameter.", tex_label = "\\gamma")
    m <= parameter(:ν, 1.0, fixed = true,
                   description = "Inverse Frisch elasticity of labor supply.", tex_label = "\\nu")
    m <= parameter(:abar, -0.5, fixed = true,
                   description = "Borrowing floor.", tex_label = "\\bar{a}")

    # More non-steadystate parameters
    m <= parameter(:ρ_tay, 0.5, (1e-3, 0.999), (1e-3, 0.999), SquareRoot(), BetaAlt(0.5, 0.2),
                   fixed = true,
                   description = "ρ_tay: Persistence in the taylor rule",
                   tex_label = "rho_{tay}")
    m <= parameter(:κ, 1.0, (1e-3, 0.999), (1e-3, 0.999), SquareRoot(), Distributions.Uniform(0.0, 1.0),
                   fixed = true,
                   description = "κ: The slope of the Phillips curve",
                   tex_label = "\\kappa")
    m <= parameter(:phipi, 1.5, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25),
                   fixed = true,
                   description = "phipi: The slope of the taylor rule",
                   tex_label = "\\phi_\\pi")

    # More steadystate parameters
    m <= parameter(:μ_s, 0., fixed = true, description = "μ_s: Mu of log normal in income")
    m <= parameter(:σ_s, 0.1, fixed = true,
                   description = "σ_s: Sigma of log normal in income")
    m <= parameter(:e_y, 1e-3, fixed = true, description = "e_y: Measurement error on GDP",
                   tex_label = "e_y")
end

"""
```
ss6!(m::RealBondMkup)
```

Initializes subspec 6 of `RealBondMkup`.
This fixes only the steady-state parameters. Used for estimating all non-steady state parameters. Right not, it's the same as ss0 since by default (when constructing the model), we do this.
"""
function ss6!(m::RealBondMkup)
    # Steady State Parameters
    m <= parameter(:R, 1.04, fixed = true,
                   description = "R: Steady-state gross real interest rate.", tex_label = "R")
    m <= parameter(:γ, 1.0, fixed = true,
                   description = "γ: CRRA Parameter.", tex_label = "\\gamma")
    m <= parameter(:ν, 1.0, fixed = true,
                   description = "Inverse Frisch elasticity of labor supply.", tex_label = "\\nu")
    m <= parameter(:abar, -0.5, fixed = true,
                   description = "Borrowing floor.", tex_label = "\\bar{a}")

    # More steadystate parameters
    m <= parameter(:μ_s, 0., fixed = true, description = "μ_s: Mu of log normal in income")
    m <= parameter(:σ_s, 0.1, fixed = true,
                   description = "σ_s: Sigma of log normal in income")
    m <= parameter(:e_y, 1e-3, fixed = true, description = "e_y: Measurement error on GDP",
                   tex_label = "e_y")
end
