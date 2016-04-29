"""
```

filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}, sys::Vector{System},
                                  z0::Vector{S}=[], vz0::Matrix{S}=[];
                                  lead::Int, Ny0::Int =0, allout::Bool = false,
                                  use_expected_rate_data::Bool = true)

filter{S<:AbstractFloat}(m::AbstractModel, data::DataFrame, sys::Vector{System},
                                  z0::Vector{S}=[], vz0::Matrix{S}=[];
                                  lead::Int, Ny0::Int =0, allout::Bool = false,
                                  use_expected_rate_data::Bool = true)

```
    
Computes and returns the filtered values of states for every state-space system in `sys`.

### Inputs

- `m`: model object
- `data`: DataFrame or matrix of data for observables
- `lead`: number of periods to forecast after the end of the data
- `sys::Vector{System}`: a vector of `System` objects specifying state-space system matrices for each draw
- `z0`: an optional `Nz x 1` initial state vector.
- `vz0`: an optional `Nz x Nz` covariance matrix of an initial state vector.
- `Ny0`: an optional scalar indicating the number of periods of
  presample (i.e. the number of periods which we don't add to the
  likelihood)
- `allout`: an optional keyword argument indicating whether we want optional output
  variables returned as well

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
                                  z0::Vector{S}=Array{S}(0),
                                  vz0::Matrix{S}=Array{S}(0,0);
                                  lead::Int=0,
                                  Ny0::Int =0,
                                  allout::Bool = false,
                                  use_expected_rate_data::Bool = true)

    
    # Convert the DataFrame to a data matrix without altering the original dataframe  
    data  = df_to_matrix(df) 
                
    filter(m, data, sys, z0, vz0, lead=lead, Ny0 =Ny0, allout=allout,
           use_expected_rate_data=use_expected_rate_data)
end


function filter{T<:AbstractModel, S<:AbstractFloat}(m::T,
                                  data::Matrix{S},
                                  sys::Vector{System{S}},   
                                  z0::Vector{S}=Array{S}(0),
                                  vz0::Matrix{S}=Array{S}(0,0);
                                  lead::Int=0,
                                  Ny0::Int =0,
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
    out   = pmap(i -> DSGE.filter(m,data,sys[i], allout=true,
                                            use_expected_rate_data=use_expected_rate_data), 1:ndraws)

    filtered_states = [x[1] for x in out]
    pred            = [x[2] for x in out]
    vpred           = [x[3] for x in out]
    
    return filtered_states, pred, vpred
end


function filter{T<:AbstractModel, S<:AbstractFloat}(m::T,
                                  data::Matrix{S},
                                  sys::System,  
                                  z0::Array{S}=Array{S}(0),
                                  vz0::Matrix{S}=Matrix{S}(0,0);
                                  lead::Int=0,
                                  Ny0::Int =0,
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

        # We have 3 regimes: presample, main sample, and expected-rate sample (starting at zlb_start_index)
        R2, R3 = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0, lead=lead, Ny0=Ny0, allout=allout)
                
        filtered_states = [R2[:filt] R3[:filt]]

        println("R2: $(size(R2[:vpred]))")
        println("R3: $(size(R3[:vpred]))")
        pred            = hcat(R2[:pred], R3[:pred])
        vpred           = cat(3, R2[:vpred], R3[:vpred])
        return filtered_states', pred, vpred
        
    else
        # regular Kalman filter with no regime-switching
        kal = kalman_filter(data', lead, CCC, TTT, DD, ZZ, VVall, z0, vz0, Ny0, allout=allout)

        return kal[:filt]', kal[:pred], kal[:vpred]
    end
end

