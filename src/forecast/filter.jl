# Immutable types s.t. we can use map on functions with kwargs
abstract FilterOutput
immutable AllOut<:FilterOutput end
immutable MinimumOut<:FilterOutput end

abstract FilterPresample
immutable IncludePresample<:FilterPresample end
immutable ExcludePresample<:FilterPresample end

"""
```
filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); cond_type::Symbol = :none, lead::Int = 0, allout::Bool = false,
    include_presample::Bool = true)

filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); lead::Int = 0, allout::Bool = false, include_presample::Bool =
    true)
```
    
Computes and returns the filtered values of states for every state-space system in `syses`.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `syses`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `z0`: an optional `Nz` x 1 initial state vector
- `vz0`: an optional `Nz` x `Nz` covariance matrix of an initial state vector
- `allout`: an optional keyword argument indicating whether we want optional
  output variables returned as well
- `include_presample`: indicates whether to include presample periods in the
  returned vector of `Kalman` objects

### Outputs

`filter` returns a vector of `Kalman` objects, which each contain the following fields:

- `logl`: value of the average log likelihood function of the SSM under assumption that
  observation noise ϵ(t) is normally distributed
- `pred`: a `Nz` x `T+lead` matrix containing one-step predicted state vectors.
- `vpred`: a `Nz` x `Nz` x `T+lead` matrix containing mean square errors of predicted
  state vectors.
- `filt`: an optional `Nz` x `T` matrix containing filtered state vectors.
- `vfilt`: an optional `Nz` x `Nz` x `T` matrix containing mean square errors of filtered
  state vectors.
"""
function filter{S<:AbstractFloat}(m::AbstractModel, df::DataFrame, syses::Vector{System{S}},
                                  z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                  cond_type::Symbol = :none, lead::Int = 0, allout::Bool = false,
                                  include_presample::Bool = true)
    
    # Convert the DataFrame to a data matrix without altering the original dataframe  
    data = df_to_matrix(m, df; cond_type = cond_type)
    filter(m, data, syses, z0, vz0; lead = lead, allout = allout, include_presample = include_presample)
end

function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, syses::Vector{System{S}},
                                  z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                  lead::Int = 0, allout::Bool = false, include_presample::Bool = true)

    # numbers of useful things
    ndraws = size(syses, 1)

    # Broadcast models and data matrices 
    models = fill(m, ndraws)
    datas = fill(data, ndraws)
    z0s = fill(z0, ndraws)
    vz0s = fill(vz0, ndraws)
    allouts = if allout
        fill(AllOut(), ndraws)
    else
        fill(MinimumOut(), ndraws)
    end
    include_presamples = if include_presample
        fill(IncludePresample(), ndraws)
    else
        fill(ExcludePresample(), ndraws)
    end
    
    # Call filter over all draws
    if use_parallel_workers(m) && nworkers() > 1
        mapfcn = pmap
    else
        mapfcn = map
    end    

    kals = mapfcn(DSGE.tricky_filter, allouts, include_presamples, models, datas, syses, z0s, vz0s)
    
    return [kal::Kalman{S} for kal in kals]
end

tricky_filter(::AllOut, ::IncludePresample, m::AbstractModel, data::Matrix, sys::System, z0::Vector, vz0::Matrix) =
    filter(m, data, sys, z0, vz0; allout = true, include_presample = true)
tricky_filter(::AllOut, ::ExcludePresample, m::AbstractModel, data::Matrix, sys::System, z0::Vector, vz0::Matrix) =
    filter(m, data, sys, z0, vz0; allout = true, include_presample = false)
tricky_filter(::MinimumOut, ::IncludePresample, m::AbstractModel, data::Matrix, sys::System, z0::Vector, vz0::Matrix) = 
    filter(m, data, sys, z0, vz0; allout = false, include_presample = true)
tricky_filter(::MinimumOut, ::ExcludePresample, m::AbstractModel, data::Matrix, sys::System, z0::Vector, vz0::Matrix) = 
    filter(m, data, sys, z0, vz0; allout = false, include_presample = false)
    
function filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::System{S},
                                  z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                  lead::Int = 0, allout::Bool = false, include_presample::Bool = true)
    
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
        kal, _, _, _ = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0;
            lead = lead, allout = allout, include_presample = include_presample)
    else
        # regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, TTT, CCC, ZZ, DD, VVall, z0, vz0;
            lead = lead, allout = allout, include_presample = include_presample)
    end

    return kal
end

"""
```
filterandsmooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); lead::Int = 0, allout::Bool = false, include_presample::Bool =
    true)

filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    syses::Vector{System{S}}, z0::Vector{S} = Vector{S}(), vz0::Matrix{S} =
    Matrix{S}(); lead::Int = 0, allout::Bool = false, include_presample::Bool =
    true)
```
    
Computes and returns the smoothed states, shocks, and pseudo-observables, as
well as the Kalman filter outputs, for every state-space system in `syses`.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `syses`: a vector of `System` objects specifying state-space
  system matrices for each draw
- `z0`: an optional `Nz` x 1 initial state vector
- `vz0`: an optional `Nz` x `Nz` covariance matrix of an initial state vector

### Outputs

- `states`: 3-dimensional array of size `nstates` x `hist_periods` x `ndraws`
  consisting of smoothed states for each draw
- `shocks`: 3-dimensional array of size `nshocks` x `hist_periods` x `ndraws`
  consisting of smoothed shocks for each draw
- `pseudo`: 3-dimensional array of size `npseudo` x `hist_periods` x `ndraws`
  consisting of pseudo-observables computed from the smoothed states for each
  draw
- `kals`: vector of Kalman objects, of length `ndraws`

where `states` and `shocks` are returned from the smoother specified by
`smoother_flag(m)`.
"""
function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
                                           syses::Vector{System{S}},
                                           z0::Vector{S} = Vector{S}(),
                                           vz0::Matrix{S} = Matrix{S}();
                                           cond_type::Symbol = :none,
                                           lead::Int = 0, allout::Bool = false)

    data = df_to_matrix(m, df; cond_type = cond_type)
    filterandsmooth(m, data, syses, z0, vz0; lead = lead)
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
                                           syses::Vector{System{S}},
                                           z0::Vector{S} = Vector{S}(),
                                           vz0::Matrix{S} = Matrix{S}();
                                           lead::Int = 0)
    # numbers of useful things
    ndraws = length(syses)

    # Broadcast models and data matrices 
    models = fill(m, ndraws)
    datas = fill(data, ndraws)
    z0s = fill(z0, ndraws)
    vz0s = fill(vz0, ndraws)
    
    # Call filter over all draws
    if use_parallel_workers(m) && nworkers() > 1
        mapfcn = pmap
    else
        mapfcn = map
    end    
    out = mapfcn(filterandsmooth, models, datas, syses, z0s, vz0s)

    # Unpack returned vector of tuples
    states = [x[1]::Matrix{S} for x in out]
    shocks = [x[2]::Matrix{S} for x in out]
    pseudo = [x[3]::Matrix{S} for x in out]
    kals   = [x[4]::Kalman{S} for x in out]

    # Splat vectors of matrices into 3-D arrays
    states = cat(3, states...)
    shocks = cat(3, shocks...)
    pseudo = cat(3, pseudo...)

    return states, shocks, pseudo, kals
end

function filterandsmooth{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::System{S},
                                           z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}();
                                           lead::Int = 0)
    ## 1. Filter

    # pull out the elements of sys
    TTT   = sys[:TTT]
    RRR   = sys[:RRR]
    CCC   = sys[:CCC]
    QQ    = sys[:QQ]
    ZZ    = sys[:ZZ]
    DD    = sys[:DD]
    VVall = sys[:VVall]

    # Call the appropriate version of the Kalman filter
    if n_anticipated_shocks(m) > 0

        # We have 3 regimes: presample, main sample, and expected-rate sample
        # (starting at index_zlb_start)
        kal, _, _, _ = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0; lead =
            lead, allout = true, include_presample = true)
    else
        # regular Kalman filter with no regime-switching
        kal = kalman_filter(m, data, TTT, CCC, ZZ, DD, VVall, z0, vz0;
            lead = lead, allout = true, include_presample = true)
    end

    ## 2. Smooth

    states, shocks = if forecast_smoother(m) == :kalman
        kalman_smoother(m, data, sys, kal[:z0], kal[:vz0], kal[:pred], kal[:vpred])
    elseif forecast_smoother(m) == :durbin_koopman
        durbin_koopman_smoother(m, data, sys, kal[:z0], kal[:vz0])
    end

    ## 3. Map smoothed states to pseudo-observables
    pseudo = if forecast_pseudoobservables(m)
        
        _, pseudo_mapping = pseudo_measurement(m)
        Z_pseudo = pseudo_mapping.ZZ
        D_pseudo = pseudo_mapping.DD
        
        D_pseudo .+ Z_pseudo * states
    else
        Matrix()
    end
    
    return states, shocks, pseudo, kal
end
