"""
```
init_subspec!(m::DSGEVAR)
```
initializes a model subspecification of a DSGEVAR. We utilize
subspecifications to distinguish between possible VARs, such as
specifying the observables and the number of lags.
"""
function init_subspec!(m::DSGEVAR)
    if subspec(m) == "ss0"
        return # initializes an empty DSGEVAR
    elseif subspec(m) == "ss1"
        ss1!(m)
    end

    return m
end

function ss1!(m::DSGEVAR)
    observables = [:obs_hours, :obs_gdpdeflator, :laborshare_t, :NominalWageGrowth]
    lags        = 4
    update!(m; observables = observables, lags = lags)
end
