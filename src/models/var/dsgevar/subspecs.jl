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
    elseif subspec(m) == "ss10"
        ss10!(m)
    elseif subspec(m) == "ss11"
        ss11!(m)
    elseif subspec(m) == "ss12"
        ss12!(m)
    elseif subspec(m) == "ss13"
        ss13!(m)
    elseif subspec(m) == "ss111"
        ss111!(m)
    else
        error("DSGEVAR subspec $(subspec(m)) is not defined.")
    end

    return m
end

# Subspecs 1-9 are "toy" subspecs running simple bivariate VARs
function ss1!(m::DSGEVAR)
    observables = [:obs_hours, :obs_gdpdeflator]
    lags        = 4
    λ           = .5
    update!(m; observables = observables, lags = lags, λ = λ)
end

function ss10!(m::DSGEVAR)
    observables = [:obs_hours, :obs_gdpdeflator, :laborshare_t, :NominalWageGrowth]
    lags        = 4
    λ           = .5
    update!(m; observables = observables, lags = lags, λ = λ)
end

function ss11!(m::DSGEVAR)
    observables = [:obs_hours, :π_t, :laborshare_t, :NominalWageGrowth]
    lags        = 4
    λ           = .5
    update!(m; observables = observables, lags = lags, λ = λ)
end

function ss12!(m::DSGEVAR)
    observables = [:obs_hours, :π_t, :laborshare_t, :NominalWageGrowth, :Epi_t]
    lags        = 4
    λ           = .5
    update!(m; observables = observables, lags = lags, λ = λ)
end

function ss13!(m::DSGEVAR)
    observables = [:obs_spread, :obs_hours, :π_t, :laborshare_t, :NominalWageGrowth, :Epi_t]
    lags        = 4
    λ           = .5
    update!(m; observables = observables, lags = lags, λ = λ)
end

function ss111!(m::DSGEVAR)
    observables = [:obs_hours, :π_t, :laborshare_t, :NominalWageGrowth]
    lags        = 4
    λ           = 2.
    update!(m; observables = observables, lags = lags, λ = λ)
end
