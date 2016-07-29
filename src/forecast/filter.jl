"""
```
filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::Vector{System},
                                  z0::Vector{S} = [], vz0::Matrix{S} = [];
                                  lead::Int, Ny0::Int = 0, allout::Bool = false)

filter{S<:AbstractFloat}(m::AbstractModel, data::DataFrame, sys::Vector{System},
                                  z0::Vector{S} = [], vz0::Matrix{S} = [];
                                  lead::Int, Ny0::Int = 0, allout::Bool = false)
```
    
Computes and returns the filtered values of states for every state-space system in `sys`.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
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
function filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, sys::Vector{System{S}},
                                  z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                  lead::Int = 0, Ny0::Int = 0, allout::Bool = false)
    
    # Convert the DataFrame to a data matrix without altering the original dataframe  
    data  = df_to_matrix(m,df) 
    filter(m, data, sys, z0, vz0, lead=lead, Ny0 =Ny0, allout=allout)
end


function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::Vector{System{S}},
                                  z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                  lead::Int = 0, Ny0::Int = 0, allout::Bool = false)

    # numbers of useful things
    ndraws = size(sys, 1)

    # Call filter over all draws
    if use_parallel_workers(m) && nworkers() > 1
        mapfcn = pmap
    else
        mapfcn = map
    end
    out = mapfcn(i -> DSGE.filter(m,data,sys[i], allout=true), 1:ndraws)

    filtered_states = [Array(x[1]) for x in out]  # to make type stable
    pred            = [Array(x[2]) for x in out]  
    vpred           = [Array(x[3]) for x in out]
    
    return filtered_states, pred, vpred
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::System,
                                  z0::Vector{S}=Vector{S}(), vz0::Matrix{S}=Matrix{S}();
                                  lead::Int=0, Ny0::Int=0, allout::Bool=false)
    
    # pull out the elements of sys
    TTT    = sys[:TTT]
    RRR    = sys[:RRR]
    CCC    = sys[:CCC]
    QQ     = sys[:QQ]
    ZZ     = sys[:ZZ]
    DD     = sys[:DD]
    VVall  = sys[:VVall]

    # Call the appropriate version of the Kalman filter
    if n_anticipated_shocks(m) > 0

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        k, R1, R2, R3 = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0,
            lead = lead, Ny0 = Ny0, allout = allout, include_presample = false)
   
        return k[:filt]', k[:pred], k[:vpred], k[:zend], k[:Pend]
    else
        # regular Kalman filter with no regime-switching
        k = kalman_filter(data', lead, CCC, TTT, DD, ZZ, VVall, z0, vz0, Ny0;
            allout = allout)

        return k[:filt]', k[:pred], k[:vpred], k[:zend], k[:Pend]
    end
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
                                           sys::Vector{System{S}},
                                           z0::Vector{S}=Vector{S}(),
                                           vz0::Matrix{S}=Matrix{S}();
                                           lead::Int=0, Ny0::Int=0, allout::Bool=false)
    output = cell(size(sys)) # hack
    for (i,s) in enumerate(sys)
        output[i] = filterandsmooth(m, df, sys, z0, vz0, lead=lead, Ny0=Ny0, allout=allout)
    end
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::System,
                                           z0::Vector{S}=Vector{S}(), vz0::Matrix{S}=Matrix{S}();
                                           lead::Int=0, Ny0::Int =0, allout::Bool = false)

    ## 1. Filter

    # pull out the elements of sys
    TTT    = sys[:TTT]
    RRR    = sys[:RRR]
    CCC    = sys[:CCC]
    QQ     = sys[:QQ]
    ZZ     = sys[:ZZ]
    DD     = sys[:DD]
    VVall  = sys[:VVall]

    # Call the appropriate version of the Kalman filter
    pred, vpred = if n_anticipated_shocks(m) > 0

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        k, R1, R2, R3 = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0, lead =
            lead, Ny0 = Ny0, allout = allout, include_presample = true)

        k[:pred], k[:vpred]
    else
        # regular Kalman filter with no regime-switching
        k = kalman_filter(data', lead, CCC, TTT, DD, ZZ, VVall, z0, vz0, Ny0,
            allout = allout)

        k[:pred], k[:vpred]
    end

    ## 2. Smooth

    smoother    = smoother_flag(m)
    
    alpha_hat, eta_hat = if smoother == :kalman
        kalman_smoother(m, data, sys[:TTT], sys[:RRR], sys[:CCC],
            sys[:QQ], sys[:ZZ], sys[:DD], A0, P0, pred, vpred)
    elseif smoother == :durbin_koopman
        durbin_koopman_smoother(m, data, sys[:TTT], sys[:RRR], sys[:CCC],
            sys[:QQ], sys[:ZZ], sys[:DD], A0, P0)
    end

    return alpha_hat, eta_hat
end
