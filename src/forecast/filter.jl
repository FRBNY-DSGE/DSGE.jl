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
                                  z0::Vector{S}=[],
                                  vz0::Matrix{S}=[]; 
                                  Ny0::Int =0,
                                  allout::Bool = false,
                                  use_expected_rate_data::Bool = true)

    # numbers of useful things
    ndraws   = n_draws(m)
    @assert size(sys,1) == ndraws

    inputs = @sync @parallel (vcat) for i = 1:ndraws
        FilterInput(m, data, sys[i], z0, vz0, lead, Ny0, allout, use_expected_rate_data)
    end
    
    filtered_states = pmap(filter, inputs)

    return filtered_states
end


function filter{S<:AbstractFloat}(m::AbstractModel,
                                  data::Matrix{S},
                                  sys::System,  
                                  z0::Array{S}=[],
                                  vz0::Matrix{S}=[];
                                  lead::Int=0,
                                  Ny0::Int =0,
                                  allout::Bool = false,
                                  use_expected_rate_data = true)

    TTT    = sys[:TTT]
    RRR    = sys[:RRR]
    CCC    = sys[:CCC]
    QQ     = sys[:QQ]
    ZZ     = sys[:ZZ]
    VVall  = sys[:VVall]

    # if initial z0 is not given, set it to the 0 vector
    if isempty(z0)
        z0 = zeros(S, n_states_augmented(m) - n_anticipated_shocks(m), 1)
    end
        
    # if initial vz0 isn't given, solve for it using discrete lyapunov equation
    if isempty(vz0)
        vz0 = solve_discrete_lyapunov(TTT, RRR*QQ*RRR')
    end

    if use_expected_rate_data

        # We have 3 regimes: presample, main sample, and expected-rate sample (starting at zlb_start_index)
        R2, R3 = kalman_filter_2part(m, data, TTT, RRR, CCC, z0, vz0, lead=lead, Ny0=Ny0, allout=allout)

        filtered_states = [R2[:filt] R3[:filt]]
        return filtered_states
        
    else
        # pull out the elements of sys and call the kalman filter
        kal = kalman_filter(data, lead, CCC, TTT, DD, ZZ, VVall, z0, vz0, Ny0, allout=allout)

        return kal[:filt]
    end
end


function kalman_filter_2part{S<:AbstractFloat}(m::AbstractModel,
                                               data::Matrix{S},
                                               TTT::Matrix{S} = Matrix{S}(0,0),
                                               RRR::Matrix{S} = Matrix{S}(0,0),
                                               CCC::Matrix{S} = Matrix{S}(0,0),
                                               z0::Array{S}   = Array{S}(0),
                                               vz0::Matrix{S} = Matrix{S}(0,0);
                                               lead::Int=0,
                                               Ny0::Int =0,
                                               allout::Bool = false,
                                               mh::Bool = false,
                                               catch_errors::Bool = false)
    

    # Partition sample into three regimes, and store associated matrices:
    # - R1: presample
    # - R2: normal
    # - R3: zero lower bound and beyond
    R1 = Dict{Symbol, Matrix{S}}()
    R2 = Dict{Symbol, Matrix{S}}()
    R3 = Dict{Symbol, Matrix{S}}()
    regime_mats = [R1, R2, R3]
    #regime_likes = zeros(T, 3)

    n_T0         = n_presample_periods(m)
    n_ant        = n_anticipated_shocks(m)
    t_zlb_start  = zlb_start_index(m)
    n_obs_no_ant = n_observables(m) - n_anticipated_shocks(m)
    n_obs        = n_observables(m)
    n_exo        = n_shocks_exogenous(m)

    n_states_no_ant = n_states_augmented(m) - n_anticipated_shocks(m)
    n_states_aug    = n_states_augmented(m)
    nstates        = n_states(m)
    regime_states   = [n_states_no_ant, n_states_no_ant, n_states_aug]

    R1[:data] = data[1:n_T0, 1:n_obs_no_ant]
    R2[:data] = data[(n_T0+1):t_zlb_start-1, 1:n_obs_no_ant]
    R3[:data] = data[t_zlb_start:end, :]

    # Step 1: solution to DSGE model - delivers transition equation for the state variables
    # transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
    # If we are in Metropolis-Hastings, then any errors coming out of gensys should be caught and a -Inf
    # posterior should be returned.
    if isempty(TTT) || isempty(RRR) || isempty(CCC)
        try
            R3[:TTT], R3[:RRR], R3[:CCC] = solve(m)
        catch err
            if catch_errors && isa(err, GensysError)
                info(err.msg)
                return LIKE_NULL_OUTPUT
            else
            rethrow(err)
            end
        end
    else
        R3[:TTT], R3[:RRR], R3[:CCC] = TTT, RRR, CCC
    end

    
    # Get normal, no ZLB matrices
    state_inds = [1:(nstates-n_ant); (nstates+1):n_states_aug]
    shock_inds = 1:(n_exo-n_ant)

    R2[:TTT] = R3[:TTT][state_inds, state_inds]
    R2[:RRR] = R3[:RRR][state_inds, shock_inds]
    R2[:CCC] = R3[:CCC][state_inds, :]

    ## step 2: define the measurement equation: X_t = ZZ*S_t + D + u_t
    ## where u_t = eta_t+MM* eps_t with var(eta_t) = EE
    ## where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

    # Get measurement equation matrices set up for normal and zlb periods
    measurement_R2 = measurement(m, R2[:TTT], R2[:RRR], R2[:CCC]; shocks=false)
    measurement_R3 = measurement(m, R3[:TTT], R3[:RRR], R3[:CCC]; shocks=true)
    for d in (:ZZ, :DD, :QQ, :VVall)
        R2[d] = measurement_R2[d]
        R3[d] = measurement_R3[d]
    end

    # Presample measurement & transition equation matrices are same as normal period
    for d in (:TTT, :RRR, :QQ, :ZZ, :DD, :VVall)
        R1[d] = R2[d]
    end

    ## step 3: compute log-likelihood using Kalman filter
    ##         note that kalman_filter function assumes a transition equation written as:
    ##         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t
    ##         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
    ##         and  VV2 = cov(eps2_t,u_u) = RRR*VV
    ##         define VVall as the joint variance of the two shocks VVall = var([eps2_tu_t])

    # Run Kalman filter on presample
    R1[:A0] = if isempty(z0)
        zeros(S, n_states_no_ant, 1)
    else
        z0
    end
    R1[:P0]         = solve_discrete_lyapunov(R1[:TTT], R1[:RRR]*R1[:QQ]*R1[:RRR]')
    out             = kalman_filter(R1[:data]', 1, zeros(S, regime_states[1], 1), R1[:TTT], R1[:DD], R1[:ZZ], R1[:VVall], R1[:A0], R1[:P0], allout=allout)

    R1[:like]       = Matrix{S}(1,1)
    R1[:like][1,1]  = out[:L]
    R1[:zend]       = out[:zend]
    R1[:Pend]       = out[:Pend]

    # Run Kalman filter on normal period
    zprev           = R1[:zend]
    Pprev           = R1[:Pend]
    out             = kalman_filter(R2[:data]', 1, zeros(regime_states[2], 1), R2[:TTT], R2[:DD], R2[:ZZ], R2[:VVall], zprev, Pprev, allout=allout)

        
    R2[:like]       = Matrix{S}(1,1)
    R2[:like][1,1]  = out[:L]
    R2[:zend]       = out[:zend]
    R2[:Pend]       = out[:Pend]
    R2[:filt]       = out[:filt]

    # Run Kalman filter on ZLB period
    # This section expands the number of states to accomodate extra states for the
    # anticipated policy shocks. It does so by taking the zend and Pend for the
    # state space without anticipated policy shocks, then shoves in nant
    # zeros in the middle of zend and Pend in the location of
    # the anticipated shock entries.
    before_shocks    = 1:(nstates-n_ant)
    after_shocks_old = (nstates-n_ant+1):(n_states_aug-n_ant)
    after_shocks_new = (nstates+1):n_states_aug

    zprev = [R2[:zend][before_shocks, :];
             zeros(S, n_ant, 1);
             R2[:zend][after_shocks_old, :]]

    Pprev                                     = zeros(S, n_states_aug, n_states_aug)
    Pprev[before_shocks, before_shocks]       = R2[:Pend][before_shocks, before_shocks]
    Pprev[before_shocks, after_shocks_new]    = R2[:Pend][before_shocks, after_shocks_old]
    Pprev[after_shocks_new, before_shocks]    = R2[:Pend][after_shocks_old, before_shocks]
    Pprev[after_shocks_new, after_shocks_new] = R2[:Pend][after_shocks_old, after_shocks_old]

    out             = kalman_filter(R3[:data]', 1, zeros(regime_states[3], 1), R3[:TTT], R3[:DD], R3[:ZZ], R3[:VVall], zprev, Pprev, allout=allout)

    R3[:like]       = Matrix{S}(1,1)
    R3[:like][1,1]  = out[:L]
    R3[:zend]       = out[:zend]
    R3[:Pend]       = out[:Pend]
    R3[:filt]       = out[:filt]

    ## Return outputs from both regimes
    return R2, R3 
   
end



## Package all inputs to the Kalman Filter into a single object to
## facilitate parallelization via pmap.
"""

`FilterInput{T<:AbstractModel,S<:AbstractFloat}`

Fields
—————–
 - `m`::T`
 - `data::Matrix{S}`
 - `lead::Int` 
 - `sys::System`
 - `z0::Array{S}`
 - `vz0::Matrix{S}`
 - `Ny0::Int`
 - `allout::Bool`
 - `use_expected_rate_data::Bool`

Notes
————-
  See the help for `kalman_filter` for an explanation of each field.
"""
type FilterInput{T<:AbstractModel,S<:AbstractFloat}
    m::T
    data::Matrix{S}
    lead::Int
    sys::System
    z0::Array{S}
    vz0::Matrix{S}
    Ny0::Int
    allout::Bool
    use_expected_rate_data::Bool
end

function FilterInput{T<:AbstractFloat}(m::AbstractModel,data::Matrix{T},sys::System{T})

    nstates = n_states_augmented(m) 

    FilterInput(m,data,1,sys,zeros(nstates,1),zeros(nstates,nstates),0,true,true)
end


function filter(f::FilterInput)

    return filter(f.m, f.data, f.sys, f.z0, f.vz0, lead=f.lead, Ny0=f.Ny0, allout=f.allout, use_expected_rate_data =f.use_expected_rate_data)  
end
