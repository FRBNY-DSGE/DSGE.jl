"""
`init_subspec!(m::Model990)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::Model990)
    if subspec(m) == "ss2"
        return
    elseif subspec(m) == "ss3"
        return
    elseif subspec(m) == "ss5"
        ss5!(m)
    else
        error("This subspec is not defined.")
    end
end

"""
```
ss3!(m::Model990)
```

Initialize subspec 3 of `Model990`. This is the same as subspec 2, but with
Iskander's changes (including `betabar` defined correctly with `σ_c`, not
`σ_ω_star`).
"""
function ss3!(m::Model990)
    return
end

"""
`ss5(m::Model990)`

Initializes subspecification 5 for Model990. Specifically, fixes ι_w
and ι_p to 0 (so that intermediate goods producers who do not readjust
prices and wages in a given period do not index to inflation.)
"""
function ss5!(m::Model990)

    m <= parameter(:ι_p, 0.0, fixed=true,
                   description= "ι_p: The persistence of last period's inflation in
                   the equation that describes the intertemporal
                   change in prices for intermediate goods producers
                   who cannot adjust prices. The change in prices is a
                   geometric average of steady-state inflation
                   (π_star, with weight (1-ι_p)) and last period's
                   inflation (π_{t-1})).",
                   tex_label="\\iota_p")


    m <= parameter(:ι_w,   0.0, fixed=true,
                   description="ι_w: No description available.",
                   tex_label="\\iota_w")
end