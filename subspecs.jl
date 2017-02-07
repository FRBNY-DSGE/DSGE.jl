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
        # estimate lnb_liq and lnb_safe (constants in liquidity and safety premia)

        m <= parameter(:lnb_liq, 0.47, (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.47, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_liq_*): Liquidity premium (percent annualized).",
                       tex_label="ln(b_{liq})")


        m <= parameter(:lnb_safe, 0.26,  (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.26, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_safe_*): Safety premium (percent annualized).",
                       tex_label="ln(b_{safe})")


    elseif subspec(m) == "ss3"
        # estimate lnb_liq and lnb_safe (constants in liquidity and safety premia),
        # with iid measurement error on the BAA spread measurement error in addition to the AAA.
        # This is equivalent to subspec 2 with estimated σ_BBB.

        m <= parameter(:lnb_liq, 0.47, (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.47, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_liq_*): Liquidity premium (percent annualized).",
                       tex_label="ln(b_{liq})")


        m <= parameter(:lnb_safe, 0.26,  (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.26, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_safe_*): Safety premium (percent annualized).",
                       tex_label="ln(b_{safe})")

        m <= parameter(:σ_BBB, 0.0, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                       fixed=false,
                       description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                       tex_label="\\sigma_{BBB}")



    elseif subspec(m) == "ss4"
        # ss2, with AR(1) process for AAA spread measurement error

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

    elseif subspec(m) == "ss5"

        # ss4, with AR(1) process for BBB spread measurement error in addition to AAA

        m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                       fixed=false,
                       description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                       tex_label="\\rho_{AAA}")

        m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                       fixed=false,
                       description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                       tex_label="\\rho_{BBB}")

        m <= parameter(:σ_BBB, 0.0, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                       fixed=false,
                       description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                       tex_label="\\sigma_{BBB}")

        m <= parameter(:lnb_liq, 0.47, (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.47, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_liq_*): Liquidity premium (percent annualized).",
                       tex_label="ln(b_{liq})")


        m <= parameter(:lnb_safe, 0.26,  (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.26, 0.1),
                       fixed=false, scaling = x -> (1 + x/100)^0.25,
                       description="ln(b_safe_*): Safety premium (percent annualized).",
                       tex_label="ln(b_{safe})")

    elseif subspec(m) == "ss6"

        # ss1, with iid measurement error on BBB spread

        m <= parameter(:σ_BBB, 0.0, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                       fixed=false,
                       description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                       tex_label="\\sigma_{BBB}")


    else
        error("This subspec is not defined.")
    end
end

