"""
`init_subspec!(m::SmetsWouters)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::SmetsWouters)

    if subspec(m) == "ss0"
        # Use σ_r_m/4 as SD for all anticipated shocks
        # See measurement equation
        return
    elseif subspec(m) == "ss1"
        # Normalize s.t. sum of variances of anticipated shocks equals the
        # variance of the contemporaneous shock
        # See measurement equation
        return
    else
        error("This subspec is not defined.")
    end

end

