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
        # estimate lnb_liq and lnb_safe

        m <= parameter(:lnb_liq, 0.47, (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.47, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_liq_*): Liquidity premium (percent annualized).",
                       tex_label="ln(b_{liq})")


        m <= parameter(:lnb_safe, 0.26,  (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.26, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_safe_*): Safety premium (percent annualized).",
                       tex_label="ln(b_{safe})")


    elseif subspec(m) == "ss3"
        # double the mean on lnb_liq
        return

    elseif subspec(m) == "ss4"

        # ss2, plus estimate the ρ_AAA and put measurement error on the AAA spread
        m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                       fixed=false,
                       description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                       tex_label="\\rho_{AAA}")

        m <= parameter(:lnb_liq, 0.47, (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.47, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_liq_*): Liquidity premium (percent annualized).",
                       tex_label="ln(b_{liq})")


        m <= parameter(:lnb_safe, 0.26,  (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.26, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_safe_*): Safety premium (percent annualized).",
                       tex_label="ln(b_{safe})")

    else
        error("This subspec is not defined.")
    end
end

