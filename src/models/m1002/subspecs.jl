"""
`init_subspec!(m::Model1002)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::Model1002)
    if subspec(m) == "ss2"
        return
    elseif subspec(m) == "ss8"
        return ss8!(m)
    else
        error("This subspec is not defined.")
    end
end


"""
```
ss8!(m::Model1002)
```

Initializes subspec 8 of `Model1002`. In this subspecification, γ_gdi
and δ_gdi are fixed at 1 and 0 respectively.
"""
function ss8!(m::Model1002)

    m <= parameter(:γ_gdi,      1., (-10., 10.),     (-10., -10.),     DSGE.Untransformed(), Normal(1., 2.),            fixed=true,
                   description="γ_gdi: No description available.",
                   tex_label="\\gamma_{gdi}")
    
    m <= parameter(:δ_gdi,   0., (-10., 10.), (-10., -10.),  DSGE.Untransformed(),
                   Normal(0.00, 2.),            fixed=true,
                   description="δ_gdi: No description available.",
                   tex_label="\\delta_{gdi}")

end
