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
    m <= parameter(:σ_mon, 0.0, fixed = true,
                   description="σ_mon: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{mon}")
    m <= parameter(:ρ_mkp, 0.0, fixed = true,
                   description="ρ_mkp: Persistence of the markup shock",
                   tex_label="\\rho_{mkp}")
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
    m <= parameter(:σ_z, 0.0, fixed = true,
                   description="σ_z: The standard deviation of the technology process.",
                   tex_label="\\sigma_z")
    m <= parameter(:ρ_mon, 0.0, fixed = true,
                   description="ρ_mon: Persistence of the monetary policy shock",
                   tex_label="\\rho_{mon}")
    m <= parameter(:σ_mon, 0.0, fixed = true,
                   description="σ_mon: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{mon}")
    m <= parameter(:ρ_mkp, 0.0, fixed = true,
                   description="ρ_mkp: Persistence of the markup shock",
                   tex_label="\\rho_{mkp}")
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
    m <= parameter(:σ_z, 0.0, fixed = true,
                   description="σ_z: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_z")
    m <= parameter(:ρ_mon, 0.0, fixed = true,
                   description="ρ_mon: Persistence of the monetary policy shock",
                   tex_label="\\rho_z")
    m <= parameter(:ρ_mkp, 0.0, fixed = true,
                   description="ρ_mkp: Persistence of the markup shock",
                   tex_label="\\rho_z")
    m <= parameter(:σ_mkp, 0.0, fixed = true,
                   description="σ_mkp: The standard deviation of the markup shock.",
                   tex_label="\\sigma_{mkp}")
end
