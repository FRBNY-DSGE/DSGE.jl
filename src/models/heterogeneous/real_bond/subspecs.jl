"""
`init_subspec!(m::RealBond)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::RealBond)
    if subspec(m) == "ss0"
        return
    elseif subspec(m) == "ss1"
        return ss1!(m)
    elseif subspec(m) == "ss2"
        return ss2!(m)
    else
        error("This subspec should be a 0")
    end
end

"""
```
ss1!(m::RealBond)
```

Initializes subspec 1 of `RealBond`. In this subspecification,
The monetary shock process is shut down with σ_mon = 0., and
ρ_z is fixed so the only parameter that is free is σ_z.
"""
function ss1!(m::RealBond)
    m <= parameter(:ρ_z, 0.0, fixed = true,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
    m <= parameter(:ρmon, 0.0, fixed = true,
                   description="ρmon: AR(1) coefficient in the monetary policy process.",
                   tex_label="\\rho_z")
    m <= parameter(:σ_z, sqrt(.007), (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_z: The standard deviation of the process describing the stationary component of productivity.",
                   tex_label="\\sigma_{z}")
    m <= parameter(:σ_mon, 0.0, fixed = true,
                   description="σ_mon: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{mon}")
end

"""
```
ss2!(m::RealBond)
```

Initializes subspec 2 of `RealBond`. In this subspecification,
The technology shock process is shut down with σ_z = 0., and
ρ_z is fixed so the only parameter that is free is σ_mon.
"""
function ss2!(m::RealBond)
    m <= parameter(:ρ_z, 0.0, fixed = true,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
    m <= parameter(:ρmon, 0.0, fixed = true,
                   description="ρmon: AR(1) coefficient in the monetary policy process.",
                   tex_label="\\rho_z")
    m <= parameter(:σ_z, 0.0, fixed = true,
                   description="σ_mon: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{mon}")
    m <= parameter(:σ_mon, 0.2380, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed=false,
                   description="σ_mon: The standard deviation of the monetary policy shock.",
                   tex_label="\\sigma_{mon}")
end
