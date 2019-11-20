import Distributions: Uniform

"""
`init_subspec!(m::KrusellSmith)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::KrusellSmith)

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
ss2!(m::KrusellSmith)
```

Initializes subspec 2 of `KrusellSmith`. Change Beta distribution to be uniform.
"""
function ss2!(m::KrusellSmith)
    m <= parameter(:ρ_z, 0.95, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Uniform(0, 1), fixed=false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label="\\rho_z")
end
