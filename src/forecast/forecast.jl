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
                                    shock_distribution::Union{Distribution, Matrix{T}}=Matrix{T}(0,0),
                                    vars_to_forecast::Vector{Symbol}=[:p1, :p2])


    # numbers of useful things
    ndraws = if m.testing
        2   # TODO: make testing less hard-codey
    else
        n_draws(m)
    end
    @assert length(sys) == ndraws

    # for now, we are ignoring pseudoobservables so these can be empty
    Z_pseudo = Matrix{Float64}(12,72)
    D_pseudo = Matrix{Float64}(12,1)
    @everywhere Z_pseudo = remotecall_fetch(1, ()->Matrix{Float64}(12,72))
    @everywhere D_pseudo = remotecall_fetch(1, ()->Matrix{Float64}(12,72))

    # put the list of variables to forecast on every node
    @eval @everywhere vars=$vars_to_forecast

    # retrieve settings for forecast
    horizon  = forecast_horizons(m)
    nshocks  = n_shocks_exogenous(m)
    ndraws   = if m.testing
        2
    else
        n_draws(m)
    end

    # set up distribution of shocks if not specified
    # For now, we construct a giant vector of distirbutions of shocks and pass
    # each to computeForecast.
    #
    # TODO: refactor so that computeForecast
    # creates its own DegenerateMvNormal based on passing the QQ
    # matrix (which has already been computed/is taking up space)
    # rather than having to copy each Distribution across nodes. This will also be much more
    # space-efficient for if forecast_kill_shocks is true.

    shock_distribution = if isempty(shock_distribution)

        # Kill shocks: make a big vector of arrays of zeros
        if forecast_kill_shocks(m)
            @sync @parallel (vcat) for i = 1:ndraws
                zeros(horizon,nshocks)
            end
        else
            if forecast_tdist_shocks(m)                    ## use t-distributed shocks
                @sync @parallel (vcat) for i = 1:ndraws
                    TDist(forecast_tdist_df_val(m))
                end
            else                                           ## use normally distributed shocks
                @sync @parallel (vcat) for i = 1:ndraws
                    DegenerateMvNormal(zeros(nshocks),sqrt(sys[i][:QQ]))
                end
            end
        end

    end


    # Forecast the states
    # todo: figure out why vars isn't working
    forecasts = pmap(i -> computeForecast(sys[i][:TTT], sys[i][:RRR], sys[i][:CCC], sys[i][:ZZ],
                                          sys[i][:DD], Z_pseudo, D_pseudo,
                                          horizon, vars_to_forecast, shock_distribution[i],
                                          vec(initial_state_draws[i,:])), 1:ndraws)

    # unpack the giant vector of dictionaries that gets returned
    states      = [x["states"] for x in forecasts]
    observables = [x["observables"] for x in forecasts]
    # pseudo      = [x["pseudo_observables"] for x in forecasts]  # not yet implemented
    # shocks      = [x["shocks"] for x in forecasts]  # not yet implemented

    return  states, observables #, shocks
end



# I'm imagining that a Forecast object could be returned from
# computeForecast rather than a dictionary. It could look something like the
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
