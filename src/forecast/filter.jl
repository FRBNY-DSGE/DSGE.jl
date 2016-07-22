"""
```
filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::Vector{System},
                                  z0::Vector{S} = [], vz0::Matrix{S} = [];
                                  lead::Int, Ny0::Int = 0, allout::Bool = false,
                                  use_expected_rate_data::Bool = true)

filter{S<:AbstractFloat}(m::AbstractModel, data::DataFrame, sys::Vector{System},
                                  z0::Vector{S} = [], vz0::Matrix{S} = [];
                                  lead::Int, Ny0::Int = 0, allout::Bool = false,
                                  use_expected_rate_data::Bool = true)
```
    
Computes and returns the filtered values of states for every state-space system in `sys`.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `lead`: number of periods to forecast after the end of the data
- `sys::Vector{System}`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `z0`: an optional `Nz x 1` initial state vector
- `vz0`: an optional `Nz x Nz` covariance matrix of an initial state vector
- `Ny0`: an optional scalar indicating the number of periods of presample
  (i.e. the number of periods which we don't add to the likelihood)
- `allout`: an optional keyword argument indicating whether we want optional
  output variables returned as well

### Outputs

`filter` returns a vector of `Kalman` objects, which each contain the following fields.
  - `logl`: value of the average log likelihood function of the SSM under assumption that
    observation noise ϵ(t) is normally distributed
  - `pred`: a `Nz` x `T+lead` matrix containing one-step predicted state vectors.
  - `vpred`: a `Nz` x `Nz` x `T+lead` matrix containing mean square errors of predicted
    state vectors.
  - `filt`: an optional `Nz` x `T` matrix containing filtered state vectors.
  - `vfilt`: an optional `Nz` x `Nz` x `T` matrix containing mean square errors of filtered
    state vectors.
"""
function filter{T<:AbstractModel, S<:AbstractFloat}(m::T,
                                  df::DataFrame,
                                  sys::Vector{System{S}},   
                                  z0::Vector{S} = Array{S}(0),
                                  vz0::Matrix{S} = Array{S}(0, 0);
                                  lead::Int = 0,
                                  Ny0::Int = 0,
                                  allout::Bool = false,
                                  use_expected_rate_data::Bool = true)
    
    # Convert the DataFrame to a data matrix without altering the original dataframe  
    data  = df_to_matrix(m,df) 
                
    filter(m, data, sys, z0, vz0, lead = lead, Ny0 = Ny0, allout = allout,
           use_expected_rate_data = use_expected_rate_data)
end


function filter{T<:AbstractModel, S<:AbstractFloat}(m::T,
                                  data::Matrix{S},
                                  sys::Vector{System{S}},   
                                  z0::Vector{S} = Array{S}(0),
                                  vz0::Matrix{S} = Array{S}(0, 0);
                                  lead::Int = 0,
                                  Ny0::Int = 0,
                                  allout::Bool = false,
                                  use_expected_rate_data::Bool = true)

    # numbers of useful things
    ndraws = if m.testing
        2
    else
        n_draws(m)
    end

    @assert size(sys,1) == ndraws
    
    # Make sure the model object and the data are defined on every node
    @everywhere m    = remotecall_fetch(1, ()->m)
    @everywhere data = remotecall_fetch(1, ()->data)

    # Call filter over all draws
    out = pmap(i -> DSGE.filter(m,data,sys[i], allout=true,
                                            use_expected_rate_data=use_expected_rate_data), 1:ndraws)

    filtered_states = [Array(x[1]) for x in out]  # to make type stable
    pred            = [Array(x[2]) for x in out]  
    vpred           = [Array(x[3]) for x in out]
    
    return filtered_states, pred, vpred
end

function filter{T<:AbstractModel}(m::T,
                                  df::DataFrame,
                                  sys::System,  
                                  z0::Array{Float64}=Array{Float64}(0),
                                  vz0::Matrix{Float64}=Matrix{Float64}(0,0);
                                  lead::Int=0,
                                  Ny0::Int =0,
                                  allout::Bool = false,
                                  use_expected_rate_data = true)

    data = df_to_matrix(m,df) 

    filter(m,data,sys,z0,vz0,lead=lead,Ny0=Ny0,allout=allout,
           use_expected_rate_data=use_expected_rate_data)
    
end
    
function filter{T<:AbstractModel, S<:AbstractFloat}(m::T,
                                  data::Matrix{S},
                                  sys::System,  
                                  z0::Array{S} = Array{S}(0),
                                  vz0::Matrix{S} = Matrix{S}(0,0);
                                  lead::Int = 0,
                                  Ny0::Int = 0,
                                  allout::Bool = false,
                                  use_expected_rate_data = true)

    
    # pull out the elements of sys
    TTT    = sys[:TTT]
    RRR    = sys[:RRR]
    CCC    = sys[:CCC]
    QQ     = sys[:QQ]
    ZZ     = sys[:ZZ]
    DD     = sys[:DD]
    VVall  = sys[:VVall]

    # Call the appropriate version of the Kalman filter
    if use_expected_rate_data

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at zlb_start_index)
        R2, R3, R1 = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0,
            lead = lead, Ny0 = Ny0, allout = allout, augment_states = true)
                
        filtered_states = hcat(R2[:filt], R3[:filt])
        pred            = hcat(R2[:pred], R3[:pred])
        vpred           = cat(3, R2[:vpred], R3[:vpred])
        zend            = R3[:zend]    # final state vector is in R3
        A0              = R2[:z0]      # initial state vector is the final state vector from R1
        P0              = R2[:vz0]
   
        return filtered_states', pred, vpred, zend, P0
        
    else
        # regular Kalman filter with no regime-switching
        kal = kalman_filter(data', lead, CCC, TTT, DD, ZZ, VVall, z0, vz0, Ny0, allout = allout)

        return kal[:filt]', kal[:pred], kal[:vpred], kal[:zend], kal[:P0]
    end
end


function filterandsmooth{T<:AbstractModel, S<:AbstractFloat}(m::T,
                                                             data::Matrix{S},
                                                             sys::System,  
                                                             z0::Array{S} = Array{S}(0),
                                                             vz0::Matrix{S} = Matrix{S}(0,0);
                                                             lead::Int = 0,
                                                             Ny0::Int = 0,
                                                             allout::Bool = false,
                                                             use_expected_rate_data = true)
    
    ##############################################################################
    ## 1. Filter
    ##############################################################################

    # pull out the elements of sys
    TTT    = sys[:TTT]
    RRR    = sys[:RRR]
    CCC    = sys[:CCC]
    QQ     = sys[:QQ]
    ZZ     = sys[:ZZ]
    DD     = sys[:DD]
    VVall  = sys[:VVall]

    # Call the appropriate version of the Kalman filter
    filtered_states, pred, vpred, zend, A0, P0 = if use_expected_rate_data

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at zlb_start_index)
        R2, R3, R1 = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0, lead =
            lead, Ny0 = Ny0, allout = allout, augment_states = true)

        filtered_states = hcat(R2[:filt], R3[:filt])
        pred            = hcat(R2[:pred], R3[:pred])
        vpred           = cat(3, R2[:vpred], R3[:vpred])
        zend            = R3[:zend]    # final state vector is in R3
        A0              = R2[:z0]      # initial state vector is the final state vector from R1
        P0              = R2[:vz0]

        filtered_states', pred, vpred, zend, A0, P0
        
    else
        # regular Kalman filter with no regime-switching
        kal = kalman_filter(data', lead, CCC, TTT, DD, ZZ, VVall, z0, vz0, Ny0, allout = allout)

        kal[:filt]', kal[:pred], kal[:vpred], kal[:zend], kal[:z0], kal[:vz0]
    end


    ##############################################################################
    ## 2. Smooth
    ##############################################################################

    # extract settings from model
    n_ant_shocks  = n_anticipated_shocks(m)
    n_ant_lags    = n_anticipated_lags(m)
    n_pre_periods = n_presample_periods(m)
    sim_smooth    = simulation_smoother_flag(m)
    
    # run simulation smoother (or kalman smoother)
    smoothed = if sim_smooth
        durbin_koopman_smoother(m, data, sys[:TTT], sys[:RRR], sys[:CCC],
            sys[:QQ], sys[:ZZ], sys[:DD], P0, Ny0 = n_pre_periods)
    else
        kalman_smoother(filtered_states[1], P0, data, pred, vpred, sys[:TTT],
            sys[:RRR], sys[:QQ], sys[:ZZ], sys[:DD], n_ant_shocks, n_ant_lags, Ny0 =
            n_pre_periods)
    end

    return smoothed
    # returns a KalmanSmooth object with fields "states" and "shocks"
end
