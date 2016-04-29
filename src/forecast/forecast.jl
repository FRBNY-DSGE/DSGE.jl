"""
`forecast(...)`

Computes forecasts for all draws, given data, system matrices, matrix
of shocks or distribution of shocks, model object

Inputs
------
- `m`: model object
- `data`: matrix of data for observables
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
                                    initial_state_draws::Matrix{T}; 
                                    shock_distribution::Union{Distribution, Matrix{T}}=Matrix{T}(0,0))


    # numbers of useful things
    ndraws = if m.testing
        2
    else
        n_draws(m)
    end
    @assert length(sys) == ndraws

    # for now, we are ignoring pseudoobservables so these can be empty
    @everywhere Z_pseudo = remotecall_fetch(1, ()->[])
    @everywhere D_pseudo = remotecall_fetch(1, ()->[])

    # retrieve settings for forecast
    horizon  = forecast_horizons(m)
    nstates  = n_states(m)

    # set up distribution of shocks if not specified
    if isempty(shock_distribution)
        
        shock_distribution = if forecast_kill_shocks(m)    # Set all shocks to zero
            zeros(horizon,nstates)
        else                                               # Construct Distribution object
            if forecast_tdist_shocks(m)                    ## use t-distributed shocks
                TDist(forecast_tdist_df_val(m))
            else                                           ## use normally distributed shocks
                Normal(0,sqrt(sys[:QQ]))
            end
        end
    end

        
    # Forecast the states

    forecasts = pmap(i -> computeForecast(sys[i][TTT], sys[i][RRR], sys[i][CCC], sys[i][ZZ],
                                          sys[i][DD], Z_pseudo, D_pseudo,
                                          forecast_horizons, vars_to_forecast, shockdistribution,
                                          initial_state_draws))
    

    # unpack the giant vector of dictionaries that gets returned
    states      = [x[:states] for x in forecasts]
    observables = [x[:observables] for x in forecasts]
    # shocks      = [x[:shocks] for x in forecasts]  # not yet implemented
        
    return states, observables #, shocks
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
