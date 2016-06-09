"""
`forecast(...)`

Computes forecasts for all draws, given data, system matrices, matrix
of shocks or distribution of shocks, model object

Inputs
------
- `m`: model object
- `data`: matrix of data of observables
- `sys::Vector{System}`: a vector of `System` objects specifying
  state-space system matrices for each draw
- `z_T`: state vector in final historical period
- `shockdistribution`: a `Distribution` to draw shock values from, or
  a matrix specifying the shock innovations in each period
- `observables::Vector{Symbol}`: a list of observables to forecast

Outputs
-------
- A vector of `Forecast` objects, which contain the forecasted values of
states, observables, and pseudoobserables, as well as the values of
shock innovations.
"""
function forecast{T<:AbstractFloat}(m::AbstractModel,
                                    sys::Vector{System{T}},
                                    initial_state_draws::Vector{Vector{T}};
                                    shock_distributions::Union{Distribution, Matrix{T}}=Matrix{T}(0,0),
                                    vars_to_forecast::Vector{Symbol}=[:p1, :p2])

    ndraws = length(sys)

    # for now, we are ignoring pseudo-observables so these can be empty
    Z_pseudo = Matrix{Float64}(12,72)
    D_pseudo = Matrix{Float64}(12,1)

    # retrieve settings for forecast
    horizon  = forecast_horizons(m)
    nshocks  = n_shocks_exogenous(m)

    # Unpack everything for call to map/pmap
    TTTs     = [s[:TTT] for s in sys]
    RRRs     = [s[:RRR] for s in sys]
    CCCs     = [s[:CCC] for s in sys]
    ZZs      = [s[:ZZ] for s in sys]
    DDs      = [s[:DD] for s in sys]

    # Prepare copies of these objects due to lack of parallel broadcast functionality
    ZZps     = [Z_pseudo for i in 1:ndraws]
    DDps     = [D_pseudo for i in 1:ndraws]
    horizons = [horizon for i in 1:ndraws]
    vars     = [vars_to_forecast for i in 1:ndraws]
        
    # set up distribution of shocks if not specified
    # For now, we construct a giant vector of distirbutions of shocks and pass
    # each to compute_forecast.
    #
    # TODO: refactor so that compute_forecast
    # creates its own DegenerateMvNormal based on passing the QQ
    # matrix (which has already been computed/is taking up space)
    # rather than having to copy each Distribution across nodes. This will also be much more
    # space-efficient when forecast_kill_shocks is true.

    shock_distributions = if isempty(shock_distributions)
        if forecast_kill_shocks(m)
            [zeros(horizon, nshocks) for i in 1:ndraws]
        else
            # use t-distributed shocks
            if forecast_tdist_shocks(m)
                [Distributions.TDist(forecast_tdist_df_val(m)) for i in 1:ndraws]
            # use normally distributed shocks
            else
                shock_distributions = Vector{DSGE.DegenerateMvNormal}(ndraws)
                for i = 1:ndraws
                    shock_distributions[i] = DSGE.DegenerateMvNormal(zeros(nshocks),sqrt(sys[i][:QQ]))
                end
                shock_distributions
            end
        end
    end

    # Forecast the states
    if use_parallel_workers(m)
        mapfcn = pmap
    else
        mapfcn = map
    end

    # Go to work!
    forecasts =  mapfcn(DSGE.compute_forecast, TTTs, RRRs, CCCs, ZZs, DDs, ZZps, DDps,
                              horizons, vars, shock_distributions, initial_state_draws)

    # unpack the giant vector of dictionaries that gets returned
    states      = [x[:states]::Array{T} for x in forecasts]
    observables = [x[:observables]::Array{T} for x in forecasts]
    pseudo      = [x[:pseudo_observables]::Array{T} for x in forecasts]
    # shocks      = [x["shocks"] for x in forecasts]  # not yet implemented

    return  states, observables, pseudo #, shocks
end



# I'm imagining that a Forecast object could be returned from
# compute_forecast rather than a dictionary. It could look something like the
# type outlined below for a single draw.
# Perhaps we could also add some fancy indexing to be able to index by names of states/observables/etc.
# Question becomes: where do we store the list of observables? In the
# model object? In a separate forecastSettings vector that we also
# pass to forecast?
# immutable Forecast{T<:AbstractFloat}
#     states::Matrix{T}
#     observables::Matrix{T}
#     pseudoobservables::Matrix{T}
#     shocks::Matrix{T}
# end
