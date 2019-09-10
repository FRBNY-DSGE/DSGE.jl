"""
`init_subspec!(m::PoolModel)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::PoolModel)
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
        error("This subspec is not defined")
    end
end

"""
```
ss1!(m::PoolModel)
```

Initialize subspec 1 of `PoolModel` which uses the original Matlab predictive densities.
"""
function ss1!(m::PoolModel)
end

"""
```
ss2!(m::PoolModel)
```

Initialize subspec 2 of `PoolModel` which uses the corrected `Model805` predictive density
and the original `Model904` predictive density from Matlab.
"""
function ss2!(m::PoolModel)
end

"""
```
ss3!(m::PoolModel)
```

Initialize subspec 3 of `PoolModel` which uses the original Matlab predictive densities.
This uses Prior 2 in the paper:

```
ρ ∼ Β(0.8,0.1),
μ ∼ N(0, Φ^{-1}(0.75))
σ^2 ∼ IG(2,1)
```
"""
function ss3!(m::PoolModel)
    m <= parameter(:ρ, 0.5, (1e-5, 0.999), (1e-5,0.999), SquareRoot(), Beta(0.8,0.1),
                   fixed = false,
                   description="ρ: persistence of AR processing underlying λ.",
                   tex_label="\\rho")
    m <= parameter(:μ, 0., (-1e3, 1e3), (-1e3, 1e3), Untransformed(), Normal(0, quantile(Normal(), 0.75)),
                   fixed = false,
                   description="μ: drift of AR processing underlying λ.",
                   tex_label="\\mu")
    m <= parameter(:σ, 1., (1e-5, 10.), (1e-5, 10.), Exponential(), RootInverseGamma(4, 1/sqrt(2)),
                   fixed = false,
                   description="σ: volatility of AR processing underlying λ.",
                   tex_label="\\sigma")

end

"""
```
ss4!(m::PoolModel)
```

Initialize subspec 4 of `PoolModel` which uses the corrected `Model805` predictive density.
This uses Prior 2 in the paper:

```
ρ ∼ Β(0.8,0.1),
μ ∼ N(0, Φ^{-1}(0.75))
σ^2 ∼ IG(2,1)
```
"""
function ss4!(m::PoolModel)
    ss3!(m)
end

"""
```
ss5!(m::PoolModel)
```

Initialize subspec 5 of `PoolModel` which uses the original Matlab predictive densities.
This uses Prior 3 in the paper:

```
ρ ∼ Β(0.8,0.1),
μ ∼ N(0, Φ^{-1}(0.75))
σ^2 ∼ IG(2,1)
```
"""
function ss5!(m::PoolModel)
    m <= parameter(:ρ, 0.5, (1e-5, 0.999), (1e-5,0.999), SquareRoot(), Beta(0.8,0.1),
                   fixed = false,
                   description="ρ: persistence of AR processing underlying λ.",
                   tex_label="\\rho")
    m <= parameter(:σ, 1., (1e-5, 10.), (1e-5, 10.), Exponential(), RootInverseGamma(4, 1/sqrt(2)),
                   fixed = false,
                   description="σ: volatility of AR processing underlying λ.",
                   tex_label="\\sigma")
end

"""
```
ss6!(m::PoolModel)
```

Initialize subspec 6 of `PoolModel` which uses the corrected `Model805` predictive density.
This uses Prior 3 in the paper:

```
ρ ∼ Β(0.8,0.1),
μ ∼ N(0, Φ^{-1}(0.75))
σ^2 ∼ IG(2,1)
```
"""
function ss6!(m::PoolModel)
    ss5!(m)
end
