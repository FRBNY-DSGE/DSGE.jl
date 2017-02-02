"""
`init_subspec!(m::Model1010)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::Model1010)
    if subspec(m) == "ss1"
        return

    elseif subspec(m) == "ss2"
        # estimate the ρ_AAA and put measurement error on the AAA spread
        m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2),
                       fixed=false,
                       description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                       tex_label="\\rho_{AAA}")

        m <= parameter(:σ_AAA, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.75),
                       fixed=false,
                       description="σ_AAA: Standard deviation on the AR(1) process for measurement error on the AAA spread.",
                   tex_label="\\sigma_{AAA}")

    elseif subspec(m) == "ss3"
        return
    else
        error("This subspec is not defined.")
    end
end

