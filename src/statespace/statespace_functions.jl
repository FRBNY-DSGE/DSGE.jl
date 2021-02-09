"""
```
compute_system(m; tvis::Bool = false, verbose = :high)
compute_system_helper(m; tvis::Bool = false, verbose = :high)
```

Given the current model parameters, compute the state-space system
corresponding to model `m`. Returns a `System` or `RegimeSwitchingSystem` object.
The keyword `apply_altpolicy` indicates whether the state-space system
should reflect an alternative policy, and the keyword `tvis` indicates
whether the state-space system involves time-varying information sets.
To use `tvis = true`, at least the setting `:tvis_information_set` must
exist. See `?DSGE.compute_tvis_system` for more information about
computing state-space systems with time-varying information sets.

compute_system just runs compute_system_helper (which actually computes
the system) and then in the case of imperfect but positive credibility,
adjusts the anticipated observables and pseudo-observables' measurement
equations.
"""
function compute_system(m::AbstractDSGEModel{T}; tvis::Bool = false,
                        set_regime_eqcond_info!::Union{Function, Nothing} = nothing,
                        verbose::Symbol = :high) where {T <: Real}

    # do NOT delete the second part of the Boolean for apply_altpolicy; it is needed for
    # the scenarios code, which continues to use :alternative_policy as a Setting
    # AND for implementing imperfect awareness via uncertain_altpolicy
    apply_altpolicy = haskey(get_settings(m), :regime_eqcond_info) || haskey(get_settings(m), :alternative_policy)

    # When regime-switching, we need to do some extra work
    if haskey(get_settings(m), :regime_switching) && get_setting(m, :regime_switching)
        # Grab these settings
        has_uncertain_altpolicy = haskey(get_settings(m), :uncertain_altpolicy)
        has_uncertain_temp_altpol = haskey(get_settings(m), :uncertain_temp_altpol)
        has_regime_eqcond_info = haskey(get_settings(m), :regime_eqcond_info)
        uncertain_altpolicy = has_uncertain_altpolicy && get_setting(m, :uncertain_altpolicy)
        uncertain_temp_altpol = has_uncertain_temp_altpol && get_setting(m, :uncertain_temp_altpol)
        is_gensys2 = haskey(get_settings(m), :gensys2) && get_setting(m, :gensys2)

        # Grab regime info dictionary, if one exists
        regime_eqcond_info = has_regime_eqcond_info ? get_setting(m, :regime_eqcond_info) : Dict{Int, EqcondEntry}()

        # Ensure :alternative_policies exists as a Setting
        @assert (!uncertain_altpolicy && !uncertain_temp_altpol) ||
            haskey(get_settings(m),
                   :alternative_policies) "To use uncertain_altpolicy/uncertain_temp_altpol, the setting :alternative_policies must be set"

        # # If uncertain_temp_altpol is false, want to make sure temp altpol period is treated as certain
        # # in the case that uncertain_altpolicy is true (hence imperfect awareness is one
        if has_uncertain_temp_altpol && !uncertain_temp_altpol && has_regime_eqcond_info
            for reg in keys(regime_eqcond_info)
                if regime_eqcond_info[reg].alternative_policy.key in get_setting(m, :temporary_altpolicy_names)
                    if ismissing(regime_eqcond_info[reg].weights) && haskey(get_settings(m), :alternative_policies)
                        # We can infer the correct weight vector from alternative_policies.
                        # If :alternative_policies is not a setting and weights is missing,
                        # then nothing will happen b/c either the lack of an :alternative_policies
                        # will trigger an error, or the weights won't be needed.
                        altpol_vec = zeros(length(get_setting(m, :alternative_policies)) + 1)
                        altpol_vec[1] = 1.0
                        regime_eqcond_info[reg] = EqcondEntry(regime_eqcond_info[reg].alternative_policy, altpol_vec)
                    elseif !ismissing(regime_eqcond_info[reg].weights)
                        # if the weights vector is not missing, then we can infer the correct length from it
                        altpol_vec = zeros(length(regime_eqcond_info[reg].weights))
                        altpol_vec[1] = 1.0
                        regime_eqcond_info[reg].weights = altpol_vec
                    end
                end
            end
        end

        # Same for uncertain_altpolicy (note: this step is unnecessary for compute_system_helper,
        ## but it is helpful for combining historical and alternative policies with
        ## the right weights, e.g. in the time_varying_credibility.jl test)
        if has_uncertain_altpolicy && !uncertain_altpolicy && has_regime_eqcond_info
            for reg in keys(regime_eqcond_info)
                if regime_eqcond_info[reg].alternative_policy.key == alternative_policy(m).key # TODO: update this to allow for a Vector{Symbol}
                    # Same logic as previous block for these if conditions                     # as with temporary policies
                    if ismissing(regime_eqcond_info[reg].weights) && haskey(get_settings(m), :alternative_policies)
                        altpol_vec = zeros(length(get_setting(m, :alternative_policies)))
                        altpol_vec[1] = 1.0
                        regime_eqcond_info[reg] = EqcondEntry(regime_eqcond_info[reg].alternative_policy, altpol_vec)
                    elseif !ismissing(regime_eqcond_info[reg].weights)
                        altpol_vec = zeros(length(regime_eqcond_info[reg].weights))
                        altpol_vec[1] = 1.0
                        regime_eqcond_info[reg].weights = altpol_vec
                    end
                end
            end
        end

        if uncertain_altpolicy && any(x -> isa(x, MultiPeriodAltPolicy), get_setting(m, :alternative_policies))

            # Handle this special case on its own
            return compute_multiperiod_altpolicy_system_helper(m; tvis = tvis, verbose = verbose)
        else
            system_main = compute_system_helper(m; tvis = tvis, verbose = verbose)
        end
    else
        system_main = compute_system_helper(m; tvis = tvis, verbose = verbose)
    end

    # If correcting measurement eqs for anticipated (pseudo) observables is unnecessary
    # (eg. running Taylor or no regime switching or no uncertainty in temp altpol or altpolicy or
    # either perfect or zero credibility - invariant perfect or zero credibility in the case of time varying),
    # then return system now.
    # The !apply_altpolicy check may be problematic after refactoring altpolicy.
    if !apply_altpolicy || !haskey(get_settings(m), :regime_switching) || !get_setting(m, :regime_switching) ||
        !has_regime_eqcond_info || # if regime_eqcond_info is not defined, then no alt policies occur
        (has_uncertain_temp_altpol && !uncertain_temp_altpol && has_uncertain_altpolicy && !uncertain_altpolicy) ||
        (!has_uncertain_temp_altpol && !has_uncertain_altpolicy)

        if haskey(get_settings(m), :regime_switching) && get_setting(m, :regime_switching) && !apply_altpolicy
            if has_regime_eqcond_info
                m <= Setting(:regime_eqcond_info, regime_info_copy)
            end
            m <= Setting(:gensys2, is_gensys2)
        end

        return system_main
    end

    # Turn off these settings temporarily to get historical policy
    m <= Setting(:uncertain_altpolicy, false)
    m <= Setting(:uncertain_temp_altpol, false)

    n_altpolicies = length(first(values(get_setting(m, :regime_eqcond_info))).weights)
    system_altpolicies = Vector{AbstractSystem}(undef, n_altpolicies)

    system_altpolicies[1] = compute_system_helper(m; tvis = tvis, verbose = verbose)

    # Save these elements. No need to deepcopy b/c we replace their values in the get_settings(m) dict with different instances
    orig_regime_eqcond_info = get_setting(m, :regime_eqcond_info)
    orig_altpol             = haskey(get_settings(m), :alternative_policy) ? get_setting(m, :alternative_policy) : nothing
    orig_tvis_infoset       = haskey(get_settings(m), :tvis_information_set) ? get_setting(m, :tvis_information_set) : nothing
    orig_temp_altpol_names  = haskey(get_settings(m), :temporary_altpolicy_names) ? get_setting(m, :temporary_altpolicy_names) : nothing

    # m <= Setting(:regime_eqcond_info, Dict{Int64, EqcondEntry}()) # TODO: maybe also delete this setting
    delete!(get_settings(m), :regime_eqcond_info)
    delete!(get_settings(m), :alternative_policy)   # does nothing if alternative_policy is not a key in get_settings(m)
    delete!(get_settings(m), :tvis_information_set) # does nothing if tvis_information_set is not a key in get_settings(m)
    delete!(get_settings(m), :temporary_altpolicy_names) # does nothing if temporary_altpolicy_names is not a key in get_settings(m)
    for i in 2:n_altpolicies # loop over alternative policies
        ## With uncertain_altpolicy off (so only calculating 1 alternative policy).
        new_altpol = get_setting(m, :alternative_policies)[i - 1]
        if isa(new_altpol, AltPolicy)
            # If AltPolicy, we assume that the user only wants the permanent altpolicy system
            # AND there is no parameter regime-switching that affects the TTT matrix or CCC vector.
            # If there is parameter regime-switching or other kinds of regime-switching
            # that affect the TTT matrix/CCC vector in the alternative policies which
            # people are using to form expectations, then the user needs to pass a
            # MultiPeriodAltPolicy type (and set the gensys2 field to false if
            # calling gensys2 is unnecessary)
            if new_altpol.key == :default_policy # in this case, we can save some extra time
                m <= Setting(:regime_switching, false) # turn off regime-switching
                delete!(get_settings(m), :regime_eqcond_info)
                system_altpolicies[i] = compute_system_helper(m; tvis = false, verbose = verbose)
                m <= Setting(:regime_switching, true) # turn off regime-switching
            else
                m <= Setting(:alternative_policy, new_altpol)
                m <= Setting(:regime_switching, false) # which does not require regime-switching
                system_altpolicies[i] = compute_system_helper(m; tvis = false, verbose = verbose)
                m <= Setting(:regime_switching, true) # needs to be turned back on for other policies
                delete!(get_settings(m), :alternative_policy)
            end
        elseif isa(new_altpol, MultiPeriodAltPolicy)
            m <= Setting(:regime_eqcond_info, new_altpol.regime_eqcond_info)
            m <= Setting(:gensys2, new_altpol.gensys2)
            if !isnothing(new_altpol.temporary_altpolicy_names)
                m <= Setting(:temporary_altpolicy_names, new_altpol.temporary_altpolicy_names)
            end
            if isnothing(new_altpol.infoset)
                delete!(get_settings(m), :tvis_information_set) # if :tvis_information_set not there, then nothing happens
                system_altpolicies[i] = compute_system_helper(m; tvis = false, verbose = verbose)
            else
                m <= Setting(:tvis_information_set, new_altpol.infoset)
                system_altpolicies[i] = compute_system_helper(m; tvis = true, verbose = verbose)
            end
        else
            error("Every element in get_setting(m, :alternative_policies) must be an AltPolicy or a MultiPeriodAltPolicy")
        end
    end

    # Now add original settings back
    m <= Setting(:regime_eqcond_info, orig_regime_eqcond_info)
    m <= Setting(:gensys2, is_gensys2)
    m <= Setting(:uncertain_altpolicy, uncertain_altpolicy)
    m <= Setting(:uncertain_temp_altpol, uncertain_temp_altpol)
    if !isnothing(orig_altpol)
        m <= Setting(:alternative_policy, orig_altpol)
    end
    if !isnothing(orig_tvis_infoset)
        m <= Setting(:tvis_information_set, orig_tvis_infoset)
    end
    if !isnothing(orig_temp_altpol_names)
        m <= Setting(:temporary_altpolicy_names, orig_temp_altpol_names)
    end

    # Checks if pseudo measurement is required
    type_tuple = (typeof(m), Vector{Matrix{T}}, Vector{Matrix{T}}, Vector{Vector{T}})
    has_pseudo = hasmethod(pseudo_measurement, type_tuple) ||
    hasmethod(pseudo_measurement, (typeof(m), Matrix{T}, Matrix{T}, Vector{T}))

    has_fwd_looking_obs = haskey(get_settings(m), :forward_looking_observables)
    if has_pseudo
        has_fwd_looking_pseudo = haskey(get_settings(m), :forward_looking_pseudo_observables)
    end

    # Correct the measurement equations for anticipated observables via weighted average
    for reg in sort!(collect(keys(get_setting(m, :regime_eqcond_info))))
        new_wt = regime_eqcond_info[reg].weights
        which_is_sys = [isa(sys, System) for sys in system_altpolicies]

        if has_fwd_looking_obs
            for k in get_setting(m, :forward_looking_observables)
                system_main.measurements[reg][:ZZ][m.observables[k], :] =
                    sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:ZZ][m.observables[k], :] :
                                       system_altpolicies[i].measurements[reg][:ZZ][m.observables[k], :]) for i in 1:length(new_wt)])
                system_main.measurements[reg][:DD][m.observables[k]] =
                    sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:DD][m.observables[k]] :
                                       system_altpolicies[i].measurements[reg][:DD][m.observables[k]]) for i in 1:length(new_wt)])
            end
        else
            system_main.measurements[reg][:ZZ] .=
                sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:ZZ] :
                                   system_altpolicies[i].measurements[reg][:ZZ]) for i in 1:length(new_wt)])
            system_main.measurements[reg][:DD] .=
                sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:DD] :
                                   system_altpolicies[i].measurements[reg][:DD]) for i in 1:length(new_wt)])
        end

        if has_pseudo
            if has_fwd_looking_pseudo
                for k in get_setting(m, :forward_looking_pseudo_observables)
                    system_main.pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :] =
                        sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:ZZ_pseudo][m.pseudo_observables[k], :] :
                                           system_altpolicies[i].pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :])
                             for i in 1:length(new_wt)])
                    system_main.pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]] =
                        sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:DD_pseudo][m.pseudo_observables[k]] :
                                           system_altpolicies[i].pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]])
                             for i in 1:length(new_wt)])
                end
            else
                system_main.pseudo_measurements[reg][:ZZ_pseudo] .=
                    sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:ZZ_pseudo] :
                                       system_altpolicies[i].pseudo_measurements[reg][:ZZ_pseudo])
                         for i in 1:length(new_wt)])
                system_main.pseudo_measurements[reg][:DD_pseudo] .=
                        sum([new_wt[i] .* (which_is_sys[i] ? system_altpolicies[i][:DD_pseudo] :
                                           system_altpolicies[i].pseudo_measurements[reg][:DD_pseudo])
                             for i in 1:length(new_wt)])
            end
        end
    end

    return system_main
end

function compute_system_helper(m::AbstractDSGEModel{T}; tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}

    if tvis
        @assert haskey(get_settings(m), :tvis_information_set) "The setting :tvis_information_set is not defined"
        n_tvis = haskey(get_settings(m), :tvis_regime_eqcond_info) ? length(get_setting(m, :tvis_regime_eqcond_info)) : 1
        if n_tvis == 1 && haskey(get_settings(m), :tvis_regime_eqcond_info)
            if haskey(get_settings(m), :regime_eqcond_info)
                if get_setting(m, :tvis_regime_eqcond_info)[1] != get_setting(m, :regime_eqcond_info)
                    warn_str = "The dictionary of functions in the Setting :regime_eqcond_info does not match the one specified " *
                    "by the length-one Setting :tvis_regime_eqcond_info. Replacing :regime_eqcond_info with the dictionary " *
                    "of functions contained in :tvis_regime_eqcond_info."
                    @warn warn_str
                    m <= Setting(:regime_eqcond_info, get_setting(m, :tvis_regime_eqcond_info)[1])
                end
            else
                m <= Setting(:regime_eqcond_info, get_setting(m, :tvis_regime_eqcond_info)[1])
            end
        end

        if n_tvis > 1 # case of n_tvis = 1 handled below to avoid constructing redundant TimeVaryingInformationSetSystem
            @assert haskey(get_settings(m), :tvis_select_system) "The setting :tvis_select_system is not defined"
            tvis_sys = compute_tvis_system(m; verbose = verbose)
            transition_eqns = Transition{T}[tvis_sys[select, reg, :transition] for (reg, select) in enumerate(tvis_sys[:select])]
            return RegimeSwitchingSystem(transition_eqns, tvis_sys[:measurements], tvis_sys[:pseudo_measurements])
        end
    end

    solution_method  = get_setting(m, :solution_method)
    regime_switching = haskey(get_settings(m), :regime_switching) && get_setting(m, :regime_switching)

    # Solve model
    if regime_switching
        n_regimes      = (regime_switching && haskey(get_settings(m), :n_regimes)) ? get_setting(m, :n_regimes) : 1
        n_hist_regimes = (regime_switching && haskey(get_settings(m), :n_hist_regimes)) ? get_setting(m, :n_hist_regimes) : 1

        if solution_method == :gensys
            # Determine which regimes should use gensys2/gensys
            gensys_regimes, gensys2_regimes = compute_gensys_gensys2_regimes(m)

            # Solve!
            TTTs, RRRs, CCCs = solve(m; regime_switching = regime_switching,
                                     regimes = collect(1:n_regimes),
                                     gensys_regimes = gensys_regimes,
                                     gensys2_regimes = gensys2_regimes,
                                     verbose = verbose)

            return RegimeSwitchingSystem(m, TTTs, RRRs, CCCs, n_regimes, tvis) # handles formation of transition and measurement equations
        else
            error("Regime switching with solution algorithms other than gensys has not been implemented.")
        end
    else
        if solution_method == :gensys

            TTT, RRR, CCC = solve(m; verbose = verbose)
            transition_equation = Transition(TTT, RRR, CCC)

            # Solve measurement equation
            measurement_equation = measurement(m, TTT, RRR, CCC)

        elseif solution_method == :klein
            # Unpacking the method from solve to hang on to TTT_jump
            if m.spec == "het_dsge"
                TTT_jump, TTT_state, eu = klein(m)
            else
                TTT_jump, TTT_state, eu = klein(m)
            end
            if eu==-1
                throw(KleinError())
            end

            TTT, RRR = klein_transition_matrices(m, TTT_state, TTT_jump)
            CCC = zeros(n_model_states(m))

            if m.spec == "real_bond_mkup"
                GDPeqn = construct_GDPeqn(m, TTT_jump)
                TTT, RRR, CCC = augment_states(m, TTT, TTT_jump, RRR, CCC, GDPeqn)
                # Measurement (needs the additional TTT_jump argument)
                measurement_equation = measurement(m, TTT, TTT_jump, RRR, CCC, GDPeqn)
            elseif m.spec == "het_dsge" || m.spec == "rep_dsge"
                TTT, RRR, CCC = augment_states(m, TTT, RRR, CCC)
                measurement_equation = measurement(m, TTT, RRR, CCC)
            else
                TTT, RRR, CCC        = augment_states(m, TTT, RRR, CCC)
                measurement_equation = measurement(m, TTT, RRR, CCC)
            end

            transition_equation = Transition(TTT, RRR, CCC)

        else
            throw("solution_method provided does not exist.")
        end

        type_tuple = (typeof(m), Matrix{T}, Matrix{T}, Vector{T})
        if hasmethod(pseudo_measurement, type_tuple)
            # Solve pseudo-measurement equation
            pseudo_measurement_equation = pseudo_measurement(m, TTT, RRR, CCC)
            return System(transition_equation, measurement_equation, pseudo_measurement_equation)
        else
            return System(transition_equation, measurement_equation)
        end
    end
end

function RegimeSwitchingSystem(m::AbstractDSGEModel{T}, TTTs::Vector{<: AbstractMatrix{T}}, RRRs::Vector{<: AbstractMatrix{T}},
                               CCCs::Vector{<: AbstractVector{T}}, n_regimes::Int, tvis::Bool) where {T <: Real}

    transition_equations = Vector{Transition{T}}(undef, n_regimes)
    for i = 1:n_regimes
        transition_equations[i] = Transition(TTTs[i], RRRs[i], CCCs[i])
    end

    # Infer which measurement and pseudo-measurement equations to use
    model_type = typeof(m)
    type_tuple = (model_type, Vector{Matrix{T}}, Vector{Matrix{T}}, Vector{Vector{T}})
    has_pseudo = hasmethod(pseudo_measurement, type_tuple) ||
        hasmethod(pseudo_measurement, (model_type, Matrix{T}, Matrix{T}, Vector{T}))
    if tvis
        if hasmethod(measurement, type_tuple)
            measurement_equations = measurement(m, TTTs, RRRs, CCCs; information_set = get_setting(m, :tvis_information_set)[reg])
        else
            measurement_equations = Vector{Measurement{T}}(undef, n_regimes)
            for reg in 1:n_regimes
                measurement_equations[reg] = measurement(m, TTTs[reg], RRRs[reg], CCCs[reg], reg = reg,
                                                         TTTs = TTTs, CCCs = CCCs,
                                                         information_set = get_setting(m, :tvis_information_set)[reg])
            end
        end

        if hasmethod(pseudo_measurement, type_tuple)
            pseudo_measurement_equations = pseudo_measurement(m, TTTs, RRRs, CCCs)
        elseif hasmethod(pseudo_measurement, (typeof(m), Matrix{T}, Matrix{T}, Vector{T}))
            pseudo_measurement_equations = Vector{PseudoMeasurement{T}}(undef, n_regimes)
            for reg in 1:n_regimes
                pseudo_measurement_equations[reg] = pseudo_measurement(m, TTTs[reg], RRRs[reg], CCCs[reg],
                                                                       reg = reg, TTTs = TTTs, CCCs = CCCs,
                                                                       information_set = get_setting(m, :tvis_information_set)[reg])
            end
        end
    else
        if hasmethod(measurement, type_tuple)
            measurement_equations = measurement(m, TTTs, RRRs, CCCs)
        else
            measurement_equations = Vector{Measurement{T}}(undef, n_regimes)
            for reg in 1:n_regimes
                measurement_equations[reg] = measurement(m, TTTs[reg], RRRs[reg], CCCs[reg], reg = reg)
            end
        end

        if hasmethod(pseudo_measurement, type_tuple)
            pseudo_measurement_equations = pseudo_measurement(m, TTTs, RRRs, CCCs)
        elseif hasmethod(pseudo_measurement, (typeof(m), Matrix{T}, Matrix{T}, Vector{T}))
            pseudo_measurement_equations = Vector{PseudoMeasurement{T}}(undef, n_regimes)
            for reg in 1:n_regimes
                pseudo_measurement_equations[reg] = pseudo_measurement(m, TTTs[reg], RRRs[reg], CCCs[reg], reg = reg)
            end
        end
    end

    if has_pseudo
        return RegimeSwitchingSystem(transition_equations, measurement_equations, pseudo_measurement_equations)
    else
        return RegimeSwitchingSystem(transition_equations, measurement_equations)
    end
end

"""
```
compute_system(m::PoolModel{T})
```

Given the current model parameters, compute the state-space system
corresponding to the PoolModel model `m`.

Outputs

```
Φ: state transition function
Ψ: likelihood function, given weights on underlying models (the states) and predictive densities
F_ϵ: structural shock distribution
F_u: likelihood function measurement error distribution
F_λ: initial distribution of λ for state transition function
"""
function compute_system(m::PoolModel{T}; tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}
    Φ, F_ϵ, F_λ = transition(m)
    Ψ, F_u = measurement(m)
    return Φ, Ψ, F_ϵ, F_u, F_λ
end

"""
```
compute_multiperiod_altpolicy_system_helper(m::AbstractDSGEModel{T}; tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}
```

calculates a state-space system with imperfect awareness and multi-period alternative policies. This calculation is performed
separately to ensure computational efficiency.

Note that the endogenous states across all policies (implemented and alternative) must be the same, and any augmentations
to the state space after calling `gensys`/`gensys2` should be done within the default `augment_states` function. The user
can typically add if-else conditions within `augment_states` to handle state space augmentations for different alternative policy.
"""
function compute_multiperiod_altpolicy_system_helper(m::AbstractDSGEModel{T}; tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}

    # Set up
    solution_method = get_setting(m, :solution_method)
    @assert solution_method == :gensys "No solution algorithm except gensys is allowed with multi-period alternative policies"

    assert_str = "To use imperfect awareness with MultiPeriodAltPolicy, " *
        "the settings :uncertain_temp_altpol and :replace_eqcond must be true, and " *
        "either regime_eqcond_info or tvis_regime_eqcond_info is a setting"
    # Note that compute_multiperiod_altpolicy system_helper will be triggered if
    # only if uncertain_altpolicy is true and at least one element in
    # get_setting(m, :alternative_policies) is a MultiPeriodAltPolicy
    @assert haskey(get_settings(m), :uncertain_temp_altpol) &&
        get_setting(m, :uncertain_temp_altpol) && haskey(get_settings(m), :replace_eqcond) &&
        get_setting(m, :replace_eqcond) && (haskey(get_settings(m), :regime_eqcond_info) ||
                                            haskey(get_settings(m), :tvis_regime_eqcond_info)) assert_str

    n_regimes = get_setting(m, :n_regimes)

    if tvis
        @assert haskey(get_settings(m), :tvis_information_set) "The setting :tvis_information_set is not defined"
        n_tvis = haskey(get_settings(m), :tvis_regime_eqcond_info) ? length(get_setting(m, :tvis_regime_eqcond_info)) : 1
        if n_tvis == 1 && haskey(get_settings(m), :tvis_regime_eqcond_info)
            if haskey(get_settings(m), :regime_eqcond_info)
                if get_setting(m, :tvis_regime_eqcond_info)[1] != get_setting(m, :regime_eqcond_info)
                    warn_str = "The dictionary of functions in the Setting :regime_eqcond_info does not match the one specified " *
                    "by the length-one Setting :tvis_regime_eqcond_info. Replacing :regime_eqcond_info with the dictionary " *
                    "of functions contained in :tvis_regime_eqcond_info."
                    @warn warn_str
                    m <= Setting(:regime_eqcond_info, get_setting(m, :tvis_regime_eqcond_info)[1])
                end
            else
                m <= Setting(:regime_eqcond_info, get_setting(m, :tvis_regime_eqcond_info)[1])
            end
        end

        if n_tvis > 1 # case of n_tvis = 1 handled below to avoid constructing redundant TimeVaryingInformationSetSystem
            assert_str =  "Multiple equilibrium conditions with time-varying information sets and multi-period " *
                "alternative policies has not been implemented"
            @assert false assert_str
        end
    end

    # First calculate the perfectly credible transition matrices for the policy actually being implemented
    system_perfect_cred_totpolicies, is_altpol = perfect_cred_multiperiod_altpolicy_systems(m; n_regimes = n_regimes, verbose = verbose)
    n_tot_policies = length(system_perfect_cred_totpolicies)

    # Checks if pseudo measurement is required
    type_tuple = (typeof(m), Vector{Matrix{T}}, Vector{Matrix{T}}, Vector{Vector{T}})
    has_pseudo = hasmethod(pseudo_measurement, type_tuple) ||
    hasmethod(pseudo_measurement, (typeof(m), Matrix{T}, Matrix{T}, Vector{T}))

    has_fwd_looking_obs = haskey(get_settings(m), :forward_looking_observables)
    if has_pseudo
        has_fwd_looking_pseudo = haskey(get_settings(m), :forward_looking_pseudo_observables)
    end

    # Form measurement and pseudo-measurement equations for imperfect awareness state space system
    gensys_regimes, gensys2_regimes = compute_gensys_gensys2_regimes(m)
    TTTs, RRRs, CCCs = solve_uncertain_multiperiod_altpolicy(m, system_perfect_cred_totpolicies, is_altpol;
                                                             gensys_regimes = gensys_regimes, gensys2_regimes = gensys2_regimes,
                                                             regimes = collect(1:n_regimes), verbose = verbose)
    system_main = RegimeSwitchingSystem(m, TTTs, RRRs, CCCs, n_regimes, tvis)

    # Correct the measurement equations for anticipated observables via weighted average
    regime_eqcond_info = get_setting(m, :regime_eqcond_info)
    which_is_system = vcat(false, is_altpol) # need to add a false to the start of is_altpol for implemented policy
    for reg in sort!(collect(keys(regime_eqcond_info)))
        new_wt = regime_eqcond_info[reg].weights

        if has_fwd_looking_obs
            for k in get_setting(m, :forward_looking_observables)
                system_main.measurements[reg][:ZZ][m.observables[k], :] =
                    sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:ZZ][m.observables[k], :] :
                                       system_perfect_cred_totpolicies[i].measurements[reg][:ZZ][m.observables[k], :]) for i in 1:length(new_wt)])
                system_main.measurements[reg][:DD][m.observables[k]] =
                    sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:DD][m.observables[k]] :
                                       system_perfect_cred_totpolicies[i].measurements[reg][:DD][m.observables[k]]) for i in 1:length(new_wt)])
            end
        else
            system_main.measurements[reg][:ZZ] .=
                sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:ZZ] :
                                   system_perfect_cred_totpolicies[i].measurements[reg][:ZZ]) for i in 1:length(new_wt)])
            system_main.measurements[reg][:DD] .=
                sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:DD] :
                                   system_perfect_cred_totpolicies[i].measurements[reg][:DD]) for i in 1:length(new_wt)])
        end

        if has_pseudo
            if has_fwd_looking_pseudo
                for k in get_setting(m, :forward_looking_pseudo_observables)
                    system_main.pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :] =
                        sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:ZZ_pseudo][m.pseudo_observables[k], :] :
                                           system_perfect_cred_totpolicies[i].pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :])
                             for i in 1:length(new_wt)])
                    system_main.pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]] =
                        sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:DD_pseudo][m.pseudo_observables[k]] :
                                           system_perfect_cred_totpolicies[i].pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]])
                             for i in 1:length(new_wt)])
                end
            else
                system_main.pseudo_measurements[reg][:ZZ_pseudo] .=
                    sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:ZZ_pseudo] :
                                       system_perfect_cred_totpolicies[i].pseudo_measurements[reg][:ZZ_pseudo])
                         for i in 1:length(new_wt)])
                system_main.pseudo_measurements[reg][:DD_pseudo] .=
                        sum([new_wt[i] .* (which_is_system[i] ? system_perfect_cred_totpolicies[i][:DD_pseudo] :
                                           system_perfect_cred_totpolicies[i].pseudo_measurements[reg][:DD_pseudo])
                             for i in 1:length(new_wt)])
            end
        end
    end

    return system_main
end

function perfect_cred_multiperiod_altpolicy_systems(m::AbstractDSGEModel{T}; n_regimes::Int = get_setting(m, :n_regimes),
                                                    verbose::Symbol = :high) where {T <: Real}

    # Set up
    m        <= Setting(:uncertain_altpolicy,   false)
    m        <= Setting(:uncertain_temp_altpol, false)
    is_altpol = [isa(x, AltPolicy) for x in get_setting(m, :alternative_policies)] # used later on
    n_totpol = length(get_setting(m, :alternative_policies)) + 1 # + 1 for the implemented policy
    sys_totpolicies = Vector{Union{System, RegimeSwitchingSystem}}(undef, n_totpol)

    # Save the following elements.
    # Note deepcopy is necessary for regime_eqcond_info b/c weights will be updated to [1., 0., ...]
    # when we compute the perfectly credible implemented policy
    orig_regime_eqcond_info = deepcopy(get_setting(m, :regime_eqcond_info))
    orig_altpol             = haskey(get_settings(m), :alternative_policy) ? get_setting(m, :alternative_policy) : nothing
    orig_tvis_infoset       = haskey(get_settings(m), :tvis_information_set) ? get_setting(m, :tvis_information_set) : nothing
    orig_temp_altpol_names  = haskey(get_settings(m), :temporary_altpolicy_names) ? get_setting(m, :temporary_altpolicy_names) : nothing
    orig_temp_altpol_len    = haskey(get_settings(m), :temporary_altpolicy_length) ? get_setting(m, :temporary_altpolicy_length) : nothing
    is_gensys2              = haskey(get_settings(m), :gensys2) && get_setting(m, :gensys2)

    ## Now calculate the perfectly credible transition matrices for all alternative policies
    sys_totpolicies[1] = compute_system(m; tvis = haskey(get_settings(m), :tvis_information_set), verbose = verbose)

    # m <= Setting(:regime_eqcond_info, Dict{Int64, EqcondEntry}()) # TODO: maybe also delete this setting
    delete!(get_settings(m), :regime_eqcond_info)
    delete!(get_settings(m), :alternative_policy)   # does nothing if alternative_policy is not a key in get_settings(m)
    delete!(get_settings(m), :tvis_information_set) # does nothing if tvis_information_set is not a key in get_settings(m)
    delete!(get_settings(m), :temporary_altpolicy_names) # does nothing if temporary_altpolicy_names is not a key in get_settings(m)

    for i in 2:n_totpol # want to populate TTTs, RRRs, and CCCs, and there are n_totpol total policies
        new_altpol = get_setting(m, :alternative_policies)[i - 1] # decrement by 1 b/c :alternative_policies only includes
        if is_altpol[i - 1]                                       # policies which aren't being implemented, as is the case for is_altpol
            # If AltPolicy, we assume that the user only wants the permanent altpolicy system
            # AND there is no parameter regime-switching that affects the TTT matrix or CCC vector.
            # If there is parameter regime-switching or other kinds of regime-switching
            # that affect the TTT matrix/CCC vector in the alternative policies which
            # people are using to form expectations, then the user needs to pass a
            # MultiPeriodAltPolicy type (and set the gensys2 field to false if
            # calling gensys2 is unnecessary)
            m <= Setting(:regime_switching, false) # turn off regime-switching
            delete!(get_settings(m), :regime_eqcond_info)
            if new_altpol.key != :default_policy # in this case, we can save some extra time
                m <= Setting(:alternative_policy, new_altpol)
            end
            sys_totpolicies[i] = compute_system(m; verbose = verbose)
            m <= Setting(:regime_switching, true) # turn off regime-switching
            if new_altpol.key != :default_policy
                delete!(get_settings(m), :alternative_policy)
            end
        else # Then it's a MultiPeriodAltPolicy
            m <= Setting(:regime_eqcond_info, new_altpol.regime_eqcond_info)
            m <= Setting(:gensys2, new_altpol.gensys2)
            if !isnothing(new_altpol.temporary_altpolicy_names)
                m <= Setting(:temporary_altpolicy_names, new_altpol.temporary_altpolicy_names)
            end
            if !isnothing(new_altpol.temporary_altpolicy_length)
                m <= Setting(:temporary_altpol_length, new_altpol.temporary_altpolicy_length)
            end
            if !isnothing(new_altpol.infoset)
                m <= Setting(:tvis_information_set, new_altpol.infoset)
                sys_totpolicies[i] = compute_system(m; tvis = true, verbose = verbose)
            else
                sys_totpolicies[i] = compute_system(m; tvis = false, verbose = verbose)
            end
        end
    end

    # Now add original settings back
    m <= Setting(:regime_eqcond_info, orig_regime_eqcond_info)
    m <= Setting(:gensys2, is_gensys2)
    m <= Setting(:uncertain_altpolicy, true)
    m <= Setting(:uncertain_temp_altpol, true)
    if !isnothing(orig_altpol)
        m <= Setting(:alternative_policy, orig_altpol)
    end
    if !isnothing(orig_tvis_infoset)
        m <= Setting(:tvis_information_set, orig_tvis_infoset)
    end
    if !isnothing(orig_temp_altpol_names)
        m <= Setting(:temporary_altpolicy_names, orig_temp_altpol_names)
    end
    if !isnothing(orig_temp_altpol_len)
        m <= Setting(:temporary_altpol_length, orig_temp_altpol_len)
    end
    m <= Setting(:uncertain_altpolicy,   true)
    m <= Setting(:uncertain_temp_altpol, true)

    # return TTTs, RRRs, CCCs, is_altpol
    return sys_totpolicies, is_altpol
end

"""
```
compute_tvis_system(m; apply_altpolicy = false, verbose = :high)
```

Given the current model parameters, compute the regime-switching
state-space system corresponding to model `m` with time-varying
information sets.

Aside from settings required for regime-switching,
additional required settings that must exist in `m` are:
- `:tvis_information_set`: Vector of `UnitRange{Int}` specifying
which regimes should be included in the information set
corresponding each regime.
- `:tvis_regime_eqcond_info`: Vector of `regime_eqcond_info`
to specify different sets of equilibrium conditions, which
generate different state space systems (usually).
- `:tvis_select_system`: Vector of `Int` to specify which
state space system to use when calculating the measurement
and pseudo measurement equations. These state space systems
correspond to different sets of equilibrium conditions (usually).
"""
function compute_tvis_system(m::AbstractDSGEModel{T}; verbose::Symbol = :high) where {T <: Real}

    # TODO: update this compute_tvis_system to compute the average over forward-looking variables
    #       when using uncertain altpol and uncertain zlb. This means that, for each
    #       vector of transition equations, we need to compute the associated average of measurement equations,
    #       which requires us to calculate the perfect credibility solutions for each vector of transition equations
    #       that are meant to have uncertainty in them.

    @assert get_setting(m, :solution_method) ==
        :gensys "Currently, the solution method must be :gensys to calculate a state-space system with time-varying information sets"

    @assert haskey(get_settings(m), :replace_eqcond) "The setting :replace_eqcond must be true to calculate a state-space system with time-varying information sets"

    # :regime_dates should have the same number of possible regimes. Any differences in eqcond
    # should be specified by tvis_regime_eqcond_info
    tvis_infoset            = get_setting(m, :tvis_information_set)
    tvis_regime_eqcond_info = get_setting(m, :tvis_regime_eqcond_info)
    tvis_select             = get_setting(m, :tvis_select_system)
    regime_switching        = get_setting(m, :regime_switching)

    n_tvis         = length(tvis_regime_eqcond_info)
    n_regimes      = regime_switching && haskey(get_settings(m), :n_regimes) ? get_setting(m, :n_regimes) : 1
    n_hist_regimes = regime_switching && haskey(get_settings(m), :n_hist_regimes) ? get_setting(m, :n_hist_regimes) : 1

    apply_altpolicy = any(.!isempty.(tvis_regime_eqcond_info)) ||
        haskey(get_settings(m), :regime_eqcond_info) || (haskey(get_settings(m), :alternative_policy) &&
        get_setting(m, :alternative_policy).key != :historical)

    # Solve model
    transitions = Vector{Vector{Transition{T}}}(undef, n_tvis)
    TTTs_vec    = Vector{Vector{Matrix{T}}}(undef, n_tvis)
    RRRs_vec    = Vector{Vector{Matrix{T}}}(undef, n_tvis)
    CCCs_vec    = Vector{Vector{Vector{T}}}(undef, n_tvis)

    for (i, regime_eqcond_info) in enumerate(tvis_regime_eqcond_info) # For each set of equilibrium conditions,
        # Update regime_eqcond_info and compute implied gensys/gensys2 regimes
        if apply_altpolicy
            m <= Setting(:regime_eqcond_info, regime_eqcond_info) # calculate the implied regime-switching system
            first_gensys2_regime = minimum(collect(keys(get_setting(m, :regime_eqcond_info))))
        else
            first_gensys2_regime = n_hist_regimes + 1
        end
        last_gensys2_regime = haskey(get_settings(m), :temporary_altpol_length) ?
            first_gensys2_regime + get_setting(m, :temporary_altpol_length) : n_regimes
        if get_setting(m, :gensys2)
            gensys_regimes = [1:first_gensys2_regime-1]
            if last_gensys2_regime != n_regimes
                append!(gensys_regimes, [last_gensys2_regime+1:n_regimes])
            end
        else
            gensys_regimes = [1:n_regimes]
        end
        gensys2_regimes = [first_gensys2_regime-1:last_gensys2_regime]

        TTTs_vec[i], RRRs_vec[i], CCCs_vec[i] = solve(m; regime_switching = regime_switching,
                                                      regimes = collect(1:n_regimes),
                                                      gensys_regimes = gensys_regimes,
                                                      gensys2_regimes = gensys2_regimes, verbose = verbose)
        transitions[i] = Vector{Transition{T}}(undef, n_regimes)
        for j in 1:n_regimes # Compute vector of Transition for constructing the TimeVaryingInformationSetSystem
            transitions[i][j] = Transition(TTTs_vec[i][j], RRRs_vec[i][j], CCCs_vec[i][j])
        end
    end

    # Infer which measurement and pseudo-measurement equations to use
    measurement_eqns = Vector{Measurement{T}}(undef,       n_regimes)
    has_pseudo       = hasmethod(pseudo_measurement, (typeof(m), Matrix{T}, Matrix{T}, Vector{T}))
    if has_pseudo # Only calculate PseudoMeasurement equation if the method exists
        pseudo_measurement_eqns = Vector{PseudoMeasurement{T}}(undef, n_regimes)
    end

    for (reg, i) in enumerate(tvis_select)
        measurement_eqns[reg] = measurement(m, TTTs_vec[i][reg], RRRs_vec[i][reg], CCCs_vec[i][reg],
                                            reg = reg, TTTs = TTTs_vec[i], CCCs = CCCs_vec[i],
                                            information_set = tvis_infoset[reg])

        if has_pseudo
            pseudo_measurement_eqns[reg] = pseudo_measurement(m, TTTs_vec[i][reg], RRRs_vec[i][reg], CCCs_vec[i][reg], reg = reg,
                                                              TTTs = TTTs_vec[i], CCCs = CCCs_vec[i],
                                                              information_set = tvis_infoset[reg])
        end
    end

    if has_pseudo
        return TimeVaryingInformationSetSystem(transitions, measurement_eqns, pseudo_measurement_eqns,
                                               tvis_infoset, tvis_select)
    else
        return TimeVaryingInformationSetSystem(transitions, measurement_eqns, tvis_infoset, tvis_select)
    end
end
