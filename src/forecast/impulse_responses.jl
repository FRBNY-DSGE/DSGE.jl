"""
```
impulse_responses(m, system)

impulse_responses(system, horizon)
```

Compute impulse responses for a single draw.

### Inputs

- `m::AbstractDSGEModel`: model object
- `system::System{S}`: state-space system matrices
- `horizon::Int`: number of periods ahead to forecast
- `flip_shocks::Bool`: Whether to compute IRFs in response to a positive shock (by default the shock magnitude is a negative 1 std. shock)

where `S<:AbstractFloat`

### Outputs

- `states::Array{S, 3}`: matrix of size `nstates` x `horizon` x `nshocks` of
  state impulse response functions
- `obs::Array{S, 3}`: matrix of size `nobs` x `horizon` x `nshocks` of
  observable impulse response functions
- `pseudo::Array{S, 3}`: matrix of size `npseudo` x `horizon` x `nshocks` of
  pseudo-observable impulse response functions
"""
function impulse_responses(m::AbstractRepModel, system::System{S};
                           flip_shocks::Bool = false) where {S<:AbstractFloat}
    horizon = impulse_response_horizons(m)
    states, obs, pseudo = impulse_responses(system, horizon, flip_shocks = flip_shocks)
    return states, obs, pseudo
end

function impulse_responses_augmented(m::Union{AbstractHetModel,RepDSGEGovDebt}, system::System{S};
                                     flip_shocks::Bool = false,
                                     use_alternate_consumption = false) where {S<:AbstractFloat}
    horizon = impulse_response_horizons(m)

    TTT_jump, TTT_state, eu  = klein(m)
    if eu == -1
        throw(KleinError())
    end

    TTT, RRR = klein_transition_matrices(m, TTT_state, TTT_jump)
    CCC      = zeros(n_model_states(m))

    TTT, RRR, CCC = #if typeof(m) != RepDSGEGovDebt{Float64}
        augment_states(m, TTT, RRR, CCC)
   #= else
        TTT, RRR, CCC
    end=#

    measurement_equation = measurement(m, TTT, RRR, CCC)
    transition_equation = Transition(TTT, RRR, CCC)
    system = System(transition_equation, measurement_equation)

    states, obs, pseudo = impulse_responses(system, horizon, flip_shocks = flip_shocks)

    if typeof(m) == RepDSGEGovDebt{Float64}
        return states, obs, pseudo
    end

    state_indices_orig = stack_indices(m.endogenous_states_original, get_setting(m, :states))
    jump_indices_orig  = stack_indices(m.endogenous_states_original, get_setting(m, :jumps))

    Qx = get_setting(m, :Qx)
    Qy = get_setting(m, :Qy)

    # In this case, the length of state_indices and jump_indices seems to give the number of
    # UNNORMALIZED states/jumps however I'm not sure if this will always be the case/if this is
    # really what we want saved into these settings...so might need to modify in future.
    states_unnormalized = Array{Float64}(undef, length(state_indices_orig),
                                         horizon, size(states, 3))
    jumps_unnormalized = Array{Float64}(undef, length(jump_indices_orig), horizon, size(states, 3))

    model_states_unnormalized = Array{Float64}(undef, length(state_indices_orig) +
                                               length(jump_indices_orig) +
                                               length(m.endogenous_states_augmented), horizon,
                                               size(states, 3))
    for i in 1:size(states, 3)
        state_inds = 1:n_backward_looking_states(m)
        states_unnormalized[:, :, i] = Matrix(Qx')*states[state_inds, :, i]
        jump_inds = n_backward_looking_states(m)+1:n_backward_looking_states(m) + n_jumps(m)
        jumps_unnormalized[:, :, i] = Matrix(Qy')*states[jump_inds, :, i]

        augmented_states = states[n_backward_looking_states(m) + n_jumps(m) + 1 : end, :, i]
        model_states_unnormalized[:, :, i] = vcat(states_unnormalized[:, :, i],
                                                  jumps_unnormalized[:, :, i],
                                                  augmented_states)
    end

    if use_alternate_consumption
        endo = m.endogenous_states_original
        c_implied = (m[:ystar]/m[:g]) - m[:xstar]
        IRFC_implied = m[:ystar]/(c_implied*m[:g])*(model_states_unnormalized[endo[:y′_t], :, :] -
                                                    model_states_unnormalized[endo[:g′_t], :, :])-
                                                    (m[:xstar]/c_implied)*
        model_states_unnormalized[endo[:I′_t], :, :]

        model_states_unnormalized[1229, :, :] = IRFC_implied

        z_consumption = model_states_unnormalized[endo[:z′_t], :, :]
        g_lag = cat(0, model_states_unnormalized[endo[:g′_t], :, :], dims = 2)[:, 1:40, :]
        IRFC_implied_lag = m[:ystar]/(c_implied*m[:g])*(model_states_unnormalized[endo[:y′_t1], :, :]  - g_lag) - (m[:xstar]/c_implied)*model_states_unnormalized[endo[:I′_t1], :, :]
        IRF_observable_implied = IRFC_implied - IRFC_implied_lag .+ z_consumption
        obs[m.observables[:obs_consumption], :, :] = IRF_observable_implied
    end

    return model_states_unnormalized, obs, pseudo
end

function impulse_responses(m::Union{AbstractHetModel,RepDSGEGovDebt}, system::System{S};
                           flip_shocks::Bool = false,
                           use_augmented_states::Bool = true,
                           use_alternate_consumption::Bool = false) where {S<:AbstractFloat}
    if use_augmented_states
        return impulse_responses_augmented(m, system, flip_shocks = flip_shocks,
                                           use_alternate_consumption = use_alternate_consumption)
    end

    horizon = impulse_response_horizons(m)
    states, obs, pseudo = impulse_responses(system, horizon, flip_shocks = flip_shocks)

    if typeof(m) == RepDSGEGovDebt
        return states, obs, pseudo
    end

    Qx = get_setting(m, :Qx)
    Qy = get_setting(m, :Qy)

    state_indices_orig = stack_indices(m.endogenous_states_original, get_setting(m, :states))
    jump_indices_orig  = stack_indices(m.endogenous_states_original, get_setting(m, :jumps))

    # In this case, the length of state_indices and jump_indices seems to give the number of
    # UNNORMALIZED states/jumps however I'm not sure if this will always be the case/if this is
    # really what we want saved into these settings...so might need to modify in future.
    states_unnormalized = Array{Float64}(undef, length(state_indices_orig),
                                         horizon, size(states, 3))
    jumps_unnormalized = Array{Float64}(undef, length(jump_indices_orig), horizon, size(states, 3))
    model_states_unnormalized = Array{Float64}(undef, length(state_indices_orig) +
                                               length(jump_indices_orig), horizon, size(states, 3))
    for i in 1:size(states, 3)
        state_inds = 1:n_backward_looking_states(m)
        states_unnormalized[:, :, i] = Matrix(Qx')*states[state_inds, :, i]
        jump_inds = n_backward_looking_states(m)+1:n_backward_looking_states(m) + n_jumps(m)
        jumps_unnormalized[:, :, i] = Matrix(Qy')*states[jump_inds, :, i]

        model_states_unnormalized[:, :, i] = vcat(states_unnormalized[:, :, i],
                                                  jumps_unnormalized[:, :, i])
    end

    if use_alternate_consumption
        endo_orig = m.endogenous_states_original
        c_implied = (m[:ystar]/m[:g]) - m[:xstar]
        IRFC_implied = m[:ystar]/(c_implied*m[:g])*(model_states_unnormalized[endo_orig[:y′_t], :, :]  - model_states_unnormalized[endo_orig[:g′_t], :, :]) - (m[:xstar]/c_implied)*model_states_unnormalized[endo_orig[:I′_t], :, :]
        model_states_unnormalized = vcat(model_states_unnormalized, IRFC_implied)

        z_consumption = model_states_unnormalized[m.endogenous_states_original[:z′_t], :, :]
        IRFC_implied_lag = cat(0, IRFC_implied, dims = 2)[:, 1:40, :]
        IRF_observable_implied = IRFC_implied - IRFC_implied_lag .+ z_consumption
        obs[m.observables[:obs_consumption], :, :] = IRF_observable_implied
    end

    return model_states_unnormalized, obs, pseudo
end

function impulse_responses(system::System{S}, horizon::Int;
                           flip_shocks::Bool = false) where {S<:AbstractFloat}
    # Setup
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 1)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    states = zeros(S, nstates, horizon, nshocks)
    obs    = zeros(S, nobs,    horizon, nshocks)
    pseudo = zeros(S, npseudo, horizon, nshocks)

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    s_0 = zeros(S, nstates)

    for i = 1:nshocks
        # Isolate single shock
        shocks = zeros(S, nshocks, horizon)
        if flip_shocks
            shocks[i, 1] = sqrt(system[:QQ][i, i]) # a positive 1 s.d. shock
        else
            shocks[i, 1] = -sqrt(system[:QQ][i, i]) # a negative 1 s.d. shock
        end
        # Iterate state space forward
        states[:, :, i], obs[:, :, i], pseudo[:, :, i], _ = forecast(system, s_0, shocks)
    end

    return states, obs, pseudo
end

# Method for specifying the subset of shocks, and the size of each shock
function impulse_responses(m::AbstractDSGEModel, system::System{S},
                           horizon::Int, shock_names::Vector{Symbol},
                           shock_values::Vector{Float64}) where S<:AbstractFloat

    # Must provide a name and value for each shock
    @assert length(shock_names) == length(shock_values)

    # Setup
    exo          = m.exogenous_shocks
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 1)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    states = zeros(S, nstates, horizon, nshocks)
    obs    = zeros(S, nobs,    horizon, nshocks)
    pseudo = zeros(S, npseudo, horizon, nshocks)

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    s_0 = zeros(S, nstates)

    for (i, shock) in enumerate(shock_names)
        # Isolate single shock
        shocks = zeros(S, nshocks, horizon)
        shocks[exo[shock], 1] = shock_values[i]

        # Iterate state space forward
        states[:, :, exo[shock]], obs[:, :, exo[shock]], pseudo[:, :, exo[shock]], _ = forecast(system, s_0, shocks)
    end

    return states, obs, pseudo
end

# Method for picking a specific shock and the size of the desired shift in the state or
# observed variable using that shock. (Omitting the feature to do pseudo-observables
# now, for simplicity)
function impulse_responses(m::AbstractDSGEModel, system::System{S},
                           horizon::Int, shock_name::Symbol,
                           var_name::Symbol, var_value::Float64) where S<:AbstractFloat
    # Setup
    var_names, var_class =
    if var_name in keys(m.endogenous_states)
        m.endogenous_states, :states
    elseif var_name in keys(m.observables)
        m.observables, :obs
    end
    exo          = m.exogenous_shocks
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 1)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    states = zeros(S, nstates, horizon, nshocks)
    obs    = zeros(S, nobs,    horizon, nshocks)
    pseudo = zeros(S, npseudo, horizon, nshocks)

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    s_0 = zeros(S, nstates)

    # Isolate single shock
    shocks = zeros(S, nshocks, horizon)
    if var_class == :states
        shocks[exo[shock_name], 1] = obtain_shock_from_desired_state_value(var_value,
                                                                      var_names[var_name],
                                                                      exo[shock_name],
                                                                      system[:RRR])
    else # == :obs
        shocks[exo[shock_name], 1] = obtain_shock_from_desired_obs_value(var_value,
                                                                    var_names[var_name],
                                                                    exo[shock_name],
                                                                    system[:ZZ],
                                                                    system[:RRR])
    end

    # Iterate state space forward
    states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)

    return states, obs, pseudo
end

### Method for computing impulse responses after permuting
### the measurement matrix ZZ, with the option
### to compute the impulse response to a Cholesky identified
### shock to observables
function impulse_responses(m::AbstractRepModel, system::System{S},
                           horizon::Int, permute_mat::Matrix,
                           shocks::Vector{Float64} =
                           Vector{Float64}(undef,0);
                           cholesky_obs_shock::Bool = false,
                           flip_shocks::Bool = false,
                           get_shocks::Bool = false) where {S<:Real}
    # Copy matrices so we can put them back at the end
    ZZ_old = copy(system[:ZZ])
    DD_old = copy(system[:DD])

    # Permute matrices
    system.measurement.ZZ = permute_mat * system[:ZZ]
    system.measurement.DD = permute_mat * system[:DD]

    # Run IRF
    out = impulse_responses(m, system, horizon, shocks;
                            flip_shocks = flip_shocks,
                            cholesky_obs_shock =
                            cholesky_obs_shock,
                            observables_order =
                            convert(Vector{Int64},
                                    permute_mat *
                                    collect(values(m.observables))),
                            get_shocks =
                            get_shocks)
    states, obs, pseudo = out[1:3]
    if get_shocks
        cholesky_shock   = permute_mat * out[4]
        structural_shock = out[5]
    end

    # Permute obs back to original order
    obs = permute_mat * obs

    # Change system back
    system.measurement.ZZ = ZZ_old
    system.measurement.DD = DD_old

    if get_shocks
        return states, obs, pseudo, cholesky_shock, structural_shock
    else
        return states, obs, pseudo
    end
end

function impulse_responses(m::AbstractRepModel,
                           system::System{S}, horizon::Int64,
                           shock_vector::Vector{S};
                           flip_shocks::Bool = false,
                           cholesky_obs_shock::Bool = false,
                           observables_order::Vector{<:Int} =
                           collect(values(m.observables)),
                           get_shocks::Bool = false,
                           use_pinv::Bool = true) where {S<:Real}
    # Set up IRFs
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    nobs    = size(system[:ZZ], 1)
    npseudo = size(system[:ZZ_pseudo], 1)

    system = DSGE.zero_system_constants(system) # Set constant system matrices to 0
    s₀ = zeros(S, nstates)

    # Back out structural shocks corresponding to Cholesky shock to observables
    if cholesky_obs_shock
        # Compute Cholesky shock
        obs_cov_mat = (system[:ZZ] * system[:RRR]) * system[:QQ] *
            (system[:ZZ] * system[:RRR])'
        cholmat = cholesky((obs_cov_mat + obs_cov_mat') ./ 2).U'

        # Filter and smooth Cholesky shock
        if use_pinv
            struct_shock_vector = pinv(system[:ZZ] * system[:RRR] * sqrt.(system[:QQ])) *
                cholmat * shock_vector
        else
            df = DataFrame(date = date_forecast_start(m))
            mobs = collect(keys(m.observables))
            for (dev_i,obs_i) in enumerate(observables_order)
                df[!,mobs[obs_i]] .= deviation[dev_i]
            end

            for (k,new_i) in zip(keys(m.observables), observables_order)
                m.observables[k] = new_i
            end
            _, struct_shock_vector, _ = smooth(m, df, system,
                                               zeros(S, nstates),
                                               zeros(S, nstates, nstates),
                                               draw_states = false,
                                               include_presample = true,
                                               in_sample = false)

            for (old_i,k) in enumerate(keys(m.observables))
                m.observables[k] = old_i
            end
        end
    elseif length(shock_vector) == nshocks
        struct_shock_vector = shock_vector
    else
        if length(shock_vector) == nobs
            error("Length of shock_vector $(length(shock_vector)) does not match" *
                  " the number of exogenous shocks $nshocks. The length matches the " *
                  "number of observables, but cholesky_obs_shock was set to false." *
                  " Did you mean otherwise?")
        else
            error("Length of shock_vector $(length(shock_vector)) does not match" *
                  " the number of exogenous shocks $nshocks.")
        end
    end

    # Run IRF
    shocks_augmented = zeros(S, nshocks, horizon)
    shocks_augmented[:,1] = flip_shocks ? struct_shock_vector : -struct_shock_vector # negative shock by default
    states, obs, pseudo = forecast(system, s₀, shocks_augmented)

    if get_shocks
        deviation = cholmat * shock_vector
        return states, obs, pseudo, deviation, struct_shock_vector
    else
        return states, obs, pseudo
    end
end

# Method for specifying impulse response to a shock that maximizes the business-cycle
# variance of a shock to a given observable, whose location is
# specified by n_obs_var.
function impulse_responses(system::System{S}, horizon::Int, frequency_band::Tuple{S,S},
                           n_obs_var::Int; flip_shocks::Bool = false,
                           get_shocks::Bool = false) where {S<:Real}
    # Setup
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    nobs    = size(system[:ZZ], 1)
    npseudo = size(system[:ZZ_pseudo], 1)
    system  = zero_system_constants(system)

    # Compute business-cycle variance maximizing shock
    increment = abs(frequency_band[1] - frequency_band[2]) / 200.
    variance = zeros(S, nshocks, nshocks)
    sd_shock = system[:RRR] * sqrt.(system[:QQ]) # one standard deviation shock
    I_state = Matrix{S}(I, nstates, nstates)
    for f = frequency_band[1]:increment:frequency_band[2]
        sumG = system[:TTT] .* exp(-im * f)
        invA = system[:ZZ] * ((I_state - sumG) \ sd_shock)
        variance += real.(invA[n_obs_var,:] * invA[n_obs_var,:]') .* increment ./
            abs(frequency_band[1] - frequency_band[2])
    end
    eigout = eigen(variance)
    q = eigout.vectors[:, argmax(eigout.values)]

    # Compute IRFs
    states = zeros(S, nstates, horizon)
    shocks_augmented      = zeros(S, nshocks, horizon)
    shocks_augmented[:,1] = flip_shocks ? -(sqrt.(system[:QQ]) * q) : (sqrt.(system[:QQ]) * q) # q -> negative shock. To stick with convention, "flip_shock" should yield a positive shock.
    states, obs, pseudo = forecast(system, zeros(S, nstates), shocks_augmented)

    if get_shocks
        return states, obs, pseudo, vec(shocks_augmented[:,1])
    else
        return states, obs, pseudo
    end
end

function impulse_responses_peg(m::AbstractRepModel, system::System{S}, horizon::Int;
                           flip_shocks::Bool = false, H::Int = 0, peg::Symbol = :all_periods, real_rate = false) where {S<:AbstractFloat}
    ## H is peg horizon

    # Setup
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 1)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    shocks = m.exogenous_shocks

    states = zeros(S, nstates, horizon) #, nshocks)
    obs    = zeros(S, nobs,    horizon) #, nshocks)
    pseudo = zeros(S, npseudo, horizon) #, nshocks)

    # Set constant system matrices to 0
    system = zero_system_constants(system)

    s_0 = zeros(S, nstates)

    ### IRFs to FFR peg for H periods
    # Set anticipated shocks to peg FFR from t+1 to t+H
    FFRpeg = -0.25/4 #/400

    PsiR1 = 0 #constant
    # ZZ matrix
    PsiR2 = zeros(nstates)
    PsiR2[m.endogenous_states[:R_t]] = 1
    if real_rate
        PsiR2[m.endogenous_states[:Eπ_t]] = -1
    end
    if H == 0
        Rht = system[:RRR][:,shocks[:rm_sh]]
    else
        Rht = system[:RRR][:,vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")])] #% Columns of Rh referring to anticipated monetary shocks
    end
    bb = zeros(H+1,1)
    MH = zeros(H+1,H+1)
    for hh = 1:H+1
        if peg == :all_periods
           # @show size(FFRpeg), size(PsiR1), size(PsiR2), size(system[:TTT]), hh, size(s_0)
            #@show size(FFRpeg - PsiR1 - PsiR2'*(system[:TTT])^hh*s_0)
            bb[hh,1] = (FFRpeg - PsiR1 - PsiR2'*(system[:TTT])^hh*s_0)
            #@show bb[hh, 1]
            MH[hh,:] = PsiR2'*(system[:TTT])^(hh-1)*Rht
        elseif peg == :some_periods
            if hh < H+1
                FFRpeg_ = 0
            else
                FFRpeg_ = FFRpeg
            end
            bb[hh,1] = (FFRpeg_ - PsiR1 - PsiR2'*(system[:TTT])^hh*s_0)[1]
            MH[hh,:] = PsiR2'*(system[:TTT])^(hh-1)*Rht
        end

    end
    monshocks = MH\bb
    #@show size(MH), bb, monshocks
    #% Simulate the model
    etpeg = zeros(nshocks,horizon)
    if H == 0
        etpeg[shocks[:rm_sh],1] = monshocks[1]
    else
        etpeg[vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")]),1] = monshocks
    end
    states[:, :], obs[:, :], pseudo[:, :] = forecast(system, s_0, etpeg)
    return states, obs, pseudo
end

# Given a desired value, x, for some state s^i_1, back out the necessary value of
# ϵ^j_1 needed s.t. s^i_1 = x.
# E.g. if I want the impulse response of output (s^i_1) to be 50 basis points
# based on a shock to productivity (ϵ^j_1), how large does that shock
# to productivity need to be?
# The answer is ϵ^j_1 = x/R_{i,j}
function obtain_shock_from_desired_state_value(state_value::Float64, state_ind::Int,
                                               shock_ind::Int, RRR::Matrix{Float64})
    return state_value/RRR[state_ind, shock_ind]
end

# Similar principle to the function for states
# The answer is ϵ^j_1 = x/{Σ_{k=1}^n Z_{i, k} R_{k, j}},
# where n is the total number of states, i is the index of the observable of interest
# y^i_1, and j is the index of the shock of interest, ϵ^j_1
function obtain_shock_from_desired_obs_value(obs_value::Float64, obs_ind::Int,
                                             shock_ind::Int, ZZ::Matrix{Float64},
                                             RRR::Matrix{Float64})
    return obs_value/dot(ZZ[obs_ind, :], RRR[:, shock_ind])
end
