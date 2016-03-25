"""
```
filter(...)
```
    
Computes and returns the filtered values of states for every parameter draw.

Inputs
------

- `m`: model object
- `data`: matrix of data for observables
- `lead`: number of periods to forecast after the end of the data
- `sys::Vector{System}`: a vector of `System` objects specifying state-space system matrices for each draw
- `z0`: an optional `Nz x 1` initial state vector.
- `vz0`: an optional `Nz x Nz` covariance matrix of an initial state vector.
- `Ny0`: an optional scalar indicating the number of periods of presample (i.e. the number  of periods which we don't add to the likelihood)
- `allout`: an optional keyword argument indicating whether we want optional output
  variables returned as well

Outputs
-------
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
function filter{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{S},
                                  lead::Int,
                                  sys::Vector{System},   
                                  z0::Vector{S}=zeros(S, n_states_augmented(m) - n_anticipated_shocks(m), 1),
                                  vz0::Matrix{S}=[]; 
                                  Ny0::Int =0,
                                  allout::Bool = false,
                                  use_expected_rate_data::Bool = true)

    # numbers of useful things
    ndraws   = n_draws(m)
    @assert size(sys,1) == ndraws

    
    # Make vector of Kalman objects to return
    filtered_states = Array{Kalman}(ndraws)
    
    # Parallelize. We want to call a version of filter that operates
    # on a single draw and then call that for each element
    for i = 1:ndraws
        filtered_states[i] = filter(m, data, lead, ..., use_expected_rate_data) 
    end
        
    return filtered_states
end


function filter{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{S},
                                  sys::System,  
                                  z0::Vector{S}=zeros(S, n_states_augmented(m) - n_anticipated_shocks(m), 1),
                                  vz0::Matrix{S}=[];
                                  lead::Int=0,
                                  Ny0::Int =0,
                                  allout::Bool = false,
                                  use_expected_rate_data = true)

    TTT = sys[:TTT]
    RRR = sys[:RRR]
    QQ  = sys[:QQ]
    
    # if initial vz0 isn't given, solve for it using discrete lyapunov equation
    if vz0 == 0
        vz0 = solve_discrete_lyapunov(TTT, RRR*QQ*RRR')
    end

    if use_expected_rate_data
        # Do stuff to calculate when regime switches that involves the model object
        # End up calling filter twice (or maybe 3 times)
        
    else
        # pull out the elements of sys and call the kalman filter
        kalman_filter(data, lead, sys[:CCC], TTT, sys[:DD], sys[:ZZ], sys[:VVall], z0, vz0, Ny0, allout)
    end
end
