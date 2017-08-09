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
        return ss2!(m)
    elseif subspec(m) == "ss3"
        return ss3!(m)
    elseif subspec(m) == "ss4"
        return ss4!(m)
    elseif subspec(m) == "ss5"
        return ss5!(m)
    elseif subspec(m) == "ss6"
        return ss6!(m)
    elseif subspec(m) == "ss7"
        return ss7!(m)
    elseif subspec(m) == "ss8"
        return ss8!(m)
    elseif subspec(m) == "ss9"
        return ss9!(m)
    elseif subspec(m) == "ss10"
        return ss10!(m)
    elseif subspec(m) == "ss11"
        return ss11!(m)
    elseif subspec(m) == "ss12"
        return ss12!(m)
    elseif subspec(m) == "ss13"
        return ss13!(m)
    elseif subspec(m) == "ss14"
        return ss14!(m)
    elseif subspec(m) == "ss15"
        return ss15!(m)
    elseif subspec(m) == "ss16"
        return ss16!(m)
    elseif subspec(m) == "ss17"
        return ss17!(m)
    elseif subspec(m) == "ss18"
        return ss18!(m)
    elseif subspec(m) == "ss19"
        return ss19!(m)
    else
        error("This subspec is not defined.")
    end
end


function ss2!(m::Model1010)
    # estimate lnb_liq and lnb_safe (constants in liquidity and safety premia)

    m <= parameter(:lnb_liq, 0.47, (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.47, 0.1),
                   fixed=false, scaling = x -> (1 + x/100)^0.25,
                   description="ln(b_liq_*): Liquidity premium (percent annualized).",
                   tex_label="ln(b_{liq})")


    m <= parameter(:lnb_safe, 0.26,  (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.26, 0.1),
                   fixed=false, scaling = x -> (1 + x/100)^0.25,
                   description="ln(b_safe_*): Safety premium (percent annualized).",
                   tex_label="ln(b_{safe})")
end

function ss3!(m::Model1010)
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

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

end


function ss4!(m::Model1010)
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
end

function ss5!(m::Model1010)

    # ss4, with AR(1) process for BBB spread measurement error in addition to AAA

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
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
end

function ss6!(m::Model1010)

    # ss1, with iid measurement error on BBB spread

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")
end

function ss7!(m::Model1010)
    # ss3, with tight prior on σ_b_liqp and σ_b_safep

    m <= parameter(:lnb_liq, 0.47, (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.47, 0.1),
                   fixed=false, scaling = x -> (1 + x/100)^0.25,
                   description="ln(b_liq_*): Liquidity premium (percent annualized).",
                   tex_label="ln(b_{liq})")


    m <= parameter(:lnb_safe, 0.26,  (1e-5, 10.),   (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.26, 0.1),
                   fixed=false, scaling = x -> (1 + x/100)^0.25,
                   description="ln(b_safe_*): Safety premium (percent annualized).",
                   tex_label="ln(b_{safe})")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")


    m <= parameter(:σ_b_liqp, 0.0269, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, 0.0269, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end

function ss8!(m::Model1010)

    # Leave lnb_liq, lnb_safe fixed. AR(1) measurement error on both spreads (like ss5).
    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")
end

function ss9!(m::Model1010)

    # Leave lnb_liq, lnb_safe fixed. AR(1) measurement error on both spreads (like ss5),
    # with prior variance centered at 0.25.

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:σ_AAA, 0.25, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.25),
                   fixed=false,
                   description="σ_AAA: Standard deviation on the AR(1) process for measurement error on the AAA spread.",
                   tex_label="\\sigma_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.25, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.25),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")
end

function ss10!(m::Model1010)

    # Leave lnb_liq, lnb_safe fixed. AR(1) measurement error on both spreads (like ss5),
    # with prior variance centered at 0.25.
    # basically ss9 mixed with ss7

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:σ_AAA, 0.25, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.25),
                   fixed=false,
                   description="σ_AAA: Standard deviation on the AR(1) process for measurement error on the AAA spread.",
                   tex_label="\\sigma_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.25, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.25),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, 0.0269, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(2., sqrt(1/400)),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, 0.0269, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(2., sqrt(1/400)),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end

function ss11!(m::Model1010)

    # Leave lnb_liq, lnb_safe fixed. AR(1) measurement error on both spreads (like ss5),
    # with prior variance centered at 0.15.
    # basically ss10 adjusted for quarterly rather than annualized

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:σ_AAA, 0.15, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.25),
                   fixed=false,
                   description="σ_AAA: Standard deviation on the AR(1) process for measurement error on the AAA spread.",
                   tex_label="\\sigma_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.15, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.25),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, 0.0269, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(2., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, 0.0269, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(2., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end

function ss12!(m::Model1010)
    # ss7 with lnb_safe and lnb_liq fixed and std dev of permanent processes centered at number adjusted for quarterly

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")


    m <= parameter(:σ_b_liqp, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end

function ss13!(m::Model1010)
    # ss12 with AR(1) on measurement error

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end

function ss14!(m::Model1010)
    # ss13, but shutting down permanent components of b shocks

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, 0., (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=true,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, 0., (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=true,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end


function ss15!(m::Model1010)
    # ss5, with fixed lnb_safe and lnb_liq

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")
end

function ss16!(m::Model1010)
    # ss13, but not adjusting for maturities in the measurement equation

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end

function ss17!(m::Model1010)
    # ss13 with :σ_b_safep fixed at 0

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, 0.0, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=true,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")
end

function ss18!(m::Model1010)
    # ss13 with ρ_z_p fixed at 0.99

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")

    m <= parameter(:ρ_z_p,     0.99, (0.0, 1.0),      (0.0, 1.0), DSGE.SquareRoot(),    BetaAlt(0.5, 0.2),
                   fixed=true, description="ρ_z_p: No description available.", tex_label="\\rho_{z^p}")
end


function ss19!(m::Model1010)
    # ss19 with a tight prior on σ_z_p

    m <= parameter(:ρ_AAA, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_AAA: AR(1) coefficient in the AAA spread process.",
                   tex_label="\\rho_{AAA}")

    m <= parameter(:ρ_BBB, 0.5, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.1),
                   fixed=false,
                   description="ρ_BBB: AR(1) coefficient in the BBB spread process.",
                   tex_label="\\rho_{BBB}")

    m <= parameter(:σ_BBB, 0.1, (1e-8, 5.),(1e-8, 5.),DSGE.Exponential(),DSGE.RootInverseGamma(2., 0.10),
                   fixed=false,
                   description="σ_BBB: Standard deviation on the AR(1) process for measurement error on the BBB spread.",
                   tex_label="\\sigma_{BBB}")

    m <= parameter(:σ_b_liqp, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(),
                   DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_liqp: Standard deviation of stationary component of liquid asset preference shifter process.",
                   tex_label="\\sigma_{b^p, liq}")

    m <= parameter(:σ_b_safep, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(),
                   DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_b_safep: Standard deviation of stationary component of safe asset preference shifter process.",
                   tex_label="\\sigma_{b^p, safe}")

    m <= parameter(:ρ_z_p, 0.99, (0.0, 1.0), (0.0, 1.0), DSGE.SquareRoot(), BetaAlt(0.5, 0.2),
                   fixed=true,
                   description="ρ_z_p: No description available.",
                   tex_label="\\rho_{z^p}")

    m <= parameter(:σ_z_p, sqrt(1/400)/4, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(),
                   DSGE.RootInverseGamma(100., sqrt(1/400)/4),
                   fixed=false,
                   description="σ_z_p: No description available.",
                   tex_label="\\sigma_{z^p}")
end
