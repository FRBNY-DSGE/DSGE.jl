"""
```
compute_system(m; apply_altpolicy = false, tvis::Bool = false, verbose = :high)
compute_system_helper(m; apply_altpolicy = false, tvis::Bool = false, verbose = :high)
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
function compute_system(m::AbstractDSGEModel{T}; apply_altpolicy::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}

    if haskey(m.settings, :regime_switching) && get_setting(m, :regime_switching)
        # Set replace_eqcond_func_dict if alternative_policy but no replace_eqcond_func_dict
        fcast_reg = haskey(m.settings, :reg_post_conditional_end) ?
            max(get_setting(m, :reg_post_conditional_end), get_setting(m, :reg_forecast_start)) :
            get_setting(m, :reg_forecast_start) # Errors if no reg_forecast_start key in settings

        if apply_altpolicy && (!haskey(m.settings, :replace_eqcond_func_dict) ||
                               !haskey(get_setting(m, :replace_eqcond_func_dict), fcast_reg)) &&
                               haskey(m.settings, :alternative_policy) && get_setting(m, :alternative_policy).key != :historical
            get_setting(m, :replace_eqcond_func_dict)[fcast_reg] = get_setting(m, :alternative_policy).replace_eqcond

            # Fixing regime_dates to account for adding fcast_reg if necessary
            cond_date = haskey(m.settings, :date_conditional_end) ?
                DSGE.iterate_quarters(get_setting(m, :date_conditional_end),1) : get_setting(m, :date_forecast_start)
            fcast_date = max(get_setting(m, :date_forecast_start), cond_date)

            if !(fcast_reg in keys(get_setting(m, :regime_dates)))
                get_setting(m, :regime_dates)[fcast_reg] = fcast_date
            elseif get_setting(m, :regime_dates)[fcast_reg] != fcast_date
                regime_dated = copy(get_setting(m, :regime_dates))

                get_setting(m, :regime_dates)[fcast_reg] = fcast_date
                for k in keys(get_setting(m, :regime_dates))
                    if k >= fcast_reg
                        get_setting(m, :regime_dates)[k+1] = regime_dated[k]
                    end
                end
            end
        end # Do for regime_dates and remove apply_altpolicy from solve.

        # Getting altpolicy weights to avoid ugliness of constantly checking which one is in use
        vary_wt, altpol_wts, altpol_wts_name = haskey(m.settings, :imperfect_awareness_varying_weights) ?
            (true, get_setting(m, :imperfect_awareness_varying_weights), :imperfect_awareness_varying_weights) :
            (haskey(m.settings, :imperfect_awareness_weights) ?
             (false, get_setting(m, :imperfect_awareness_weights), :imperfect_awareness_weights) :
             (false, [2.], :irrelevant))

        # Grab these settings
        has_uncertain_altpolicy = haskey(m.settings, :uncertain_altpolicy)
        has_uncertain_zlb = haskey(m.settings, :uncertain_zlb)
        has_replace_eqcond_func = haskey(m.settings, :replace_eqcond_func_dict)
        uncertain_altpolicy = has_uncertain_altpolicy && get_setting(m, :uncertain_altpolicy)
        uncertain_zlb = has_uncertain_zlb && get_setting(m, :uncertain_zlb)
        gensys2 = haskey(m.settings, :gensys2) && get_setting(m, :gensys2)

        # Set replace_eqcond to nothing if !apply_altpolicy
        if !apply_altpolicy && has_replace_eqcond_func
            replace_eq_copy = copy(get_setting(m, :replace_eqcond_func_dict))
            delete!(m.settings, :replace_eqcond_func_dict)
            has_replace_eqcond_func = false
            m <= Setting(:gensys2, false)
        end

        # If uncertain_zlb is false, want to make sure ZLB period is treated as certain.
        if has_uncertain_zlb && !uncertain_zlb && has_replace_eqcond_func && altpol_wts_name != :irrelevant
            if vary_wt
                for reg in keys(altpol_wts)
                    if get_setting(m, :replace_eqcond_func_dict)[reg] == DSGE.zero_rate_replace_eq_entries
                        altpol_vec = zeros(length(altpol_wts[reg]))
                        altpol_vec[1] = 1.0
                        altpol_wts[reg] = altpol_vec
                    end
                end
                ## Without time-varying credibility, uncertain_zlb and uncertain_altpolicy must be the same
            elseif haskey(m.settings, :uncertain_altpolicy) && !get_setting(m, :uncertain_altpolicy)
                @warn "If not using time varying credibility, uncertain_zlb and uncertain_altpolicy must be the same."
                if get_setting(m, :replace_eqcond_func_dict)[reg] == DSGE.zero_rate_replace_eq_entries
                    altpol_vec = zeros(length(altpol_wts[reg]))
                    altpol_vec[1] = 1.0
                    regime_eqcond_info[reg].weights = altpol_vec
                end
            end
            ## Without time-varying credibility, uncertain_zlb and uncertain_altpolicy must be the same
        elseif has_uncertain_altpolicy && has_uncertain_zlb && uncertain_zlb != uncertain_altpolicy
            error("If not using time varying credibility, uncertain_zlb and uncertain_altpolicy must be the same.")
        end
    end

    # Same for uncertain_altpolicy (note: unnecessary for compute_system_helper
    ## only helpful for combining historical and alternative policies with
    ## the right weights later on
    if has_uncertain_altpolicy && !uncertain_altpolicy &&
        has_regime_eqcond_info && altpol_wts_name != :irrelevant
        if vary_wt
            for reg in keys(regime_eqcond_info)
                if regime_eqcond_info[reg].alternative_policy == alternative_policy(m)
                    altpol_vec = zeros(length(regime_eqcond_info[reg].weights))
                    altpol_vec[1] = 1.0
                    regime_eqcond_info[reg].weights = altpol_vec
                end
            end
        end # Else accounted for in ZLB
    end

    system_main = compute_system_helper(m; apply_altpolicy = apply_altpolicy, tvis = tvis, verbose = verbose)

    # If correcting measurement eqs for anticipated (pseudo) observables is unnecessary
    # (eg. running Taylor or no regime switching or no uncertainty in ZLB or altpolicy or
    # either perfect or zero credibility - invariant perfect or zero credibility in the case of time varying),
    # then return system now.
    # The !apply_altpolicy check may be problematic after refactoring altpolicy.
    if !apply_altpolicy || !haskey(m.settings, :regime_switching) || !get_setting(m, :regime_switching) ||
        !has_regime_eqcond_info || # if replace_eqcond_func_dict is not defined, then no alt policies occur
        (has_uncertain_zlb && !uncertain_zlb && has_uncertain_altpolicy && !uncertain_altpolicy)
        ## TODO: Setting names should change once refactoring done

        if !apply_altpolicy
            m <= Setting(:replace_eqcond_func_dict, replace_eq_copy)
            m <= Setting(:gensys2, gensys2)
        end
        return system_main
    end


    # Turn off these settings temporarily to get historical policy
    m <= Setting(:regime_switching, false)
    m <= Setting(:uncertain_altpolicy, false)
    m <= Setting(:uncertain_zlb, false)

    # Save this name to be replaced later
    altpol_vec_orig = get_setting(m, altpol_wts_name)
    replace_eq_copy = copy(get_setting(m, :replace_eqcond_func_dict))
    delete!(m.settings, :replace_eqcond_func_dict)
    gensys2 = haskey(m.settings, :gensys2) && get_setting(m, :gensys2)
    m <= Setting(:gensys2, false)

    # Called taylor but should run whatever is specified as historical policy
    ## either in alternative_policies or imperfect_credibility_historical_policy setting
    #=    altpols = get_setting(m, :alternative_policies)
    orig_replace_eqcond_func_dict = get_setting(m, :replace_eqcond_func_dict)
    for altpol in altpols
    m <= Setting(:replace_eqcond_func_dict, Dict(1 => altpol.replace_eq_entries)) # one possible way to get perfect credibility
    system_altpol = compute_system_helper(m; apply_altpolicy = false, tvis = tvis, verbose = verbose) # of alternative policies
    end=#

    system_taylor = compute_system_helper(m; apply_altpolicy = false, tvis = tvis, verbose = verbose)
    m <= Setting(:replace_eqcond_func_dict, replace_eq_copy)
    m <= Setting(:gensys2, gensys2)

    # n_altpolicies = vary_wt ? (length(first(values(altpol_wts))) - 1) : (length(altpol_wts) - 1)
    # n_altpolicies = vary_wt ? length(first(values(altpol_wts))) : length(altpol_wts)
    n_altpolicies = 1
    system_altpolicies = Vector{RegimeSwitchingSystem}(undef, n_altpolicies)

    m <= Setting(:regime_switching, true) # turn back on, but still keep uncertain_altpolicy, uncertain_zlb off => perfect credibility

    # TODO: THIS BLOCK OF CODE IS INEFFICIENT. The problem is this. We want to compute the perfect credibility
    # state space system for each possible alternative policy when there is imperfect awareness/policy uncertainty.
    # This code does it by setting the credibility of different policies to 1 and looping. HOWEVER,
    # this code incurs unnecessary costs b/c we compute the transition matrices every loop for all alternative policies
    # and calculate the implied transition equations given the weights vector. What we really want to do is calculate
    # the transition equations w/out the uncertain_altpolicy thing on at all and just use the correct replace_eqcond_func_dict
    # or solve/eqcond functions.
    #
    # To avoid this problem, currently, we are hard-coding that there are only two alternative policies, the desired one
    # and the historical one, which we assume to be the transition equations when there is no regime-switching,
    # i.e. the system_taylor calculated above.
    for i in 1:n_altpolicies # loop over alternative policies, noting that we've already computed the historical policy
        ## With alternative policy weights = [1., 0.] in all forecast regimes for given altpolicy
        #  Commented out here b/c uncertain_altpolicy is already off, so no need to bother with this.
        #  In general, we don't even want this approach b/c it involves calculating all alternative policies multiple times.
        #  It is more efficient to find a way to compute each altpolicy individually ONCE.
        #=        if vary_wt
        altpol_vec = deepcopy(altpol_wts) # abuse of the name here, technically. altpol_vec is a Dictionary in this case
        for j in keys(altpol_vec)
        altpol_vec[j] = zeros(length(n_altpolicies) + 1) # weight on n_altpolicies first, then weight on historical policy
        altpol_vec[j][i] = 1.0 # perfect credibility on each alternative policy
        end
        else
        altpol_vec = zeros(n_altpolicies + 1)
        altpol_vec[i] = 1.0
        end
        m <= Setting(altpol_wts_name, altpol_vec)=#
        system_altpolicies[i] = compute_system_helper(m; apply_altpolicy = apply_altpolicy, tvis = tvis, verbose = verbose)
    end

    # Now add uncertain altpolicy and zlb back
    m <= Setting(:uncertain_altpolicy, uncertain_altpolicy)
    m <= Setting(:uncertain_zlb, uncertain_zlb)
    #m <= Setting(altpol_wts_name, altpol_vec_orig)

    # Checks if pseudo measurement is required
    type_tuple = (typeof(m), Vector{Matrix{T}}, Vector{Matrix{T}}, Vector{Vector{T}})
    has_pseudo = hasmethod(pseudo_measurement, type_tuple) ||
    hasmethod(pseudo_measurement, (typeof(m), Matrix{T}, Matrix{T}, Vector{T}))

    has_fwd_looking_obs = haskey(get_settings(m), :forward_looking_observables)
    if has_pseudo
        has_fwd_looking_pseudo = haskey(get_settings(m), :forward_looking_pseudo_observables)
    end

    ## Correct the measurement equations for anticipated observables via convex combination
    for reg in sort!(collect(keys(get_setting(m, :replace_eqcond_func_dict))))
        new_wt = vary_wt ? altpol_wts[reg] : altpol_wts

        if has_fwd_looking_obs
            for k in get_setting(m, :forward_looking_observables)
                # COMMENTED CODE IS THE DESIRED VERSION in some form or fashion. May want to write a helper function which
                # figures out whether a given system is regime-switching or not.
                #=                system_main.measurements[reg][:ZZ][m.observables[k], :] =
                sum([new_wt[i] .* system_altpolicies[i].measurements[reg][:ZZ][m.observables[k], :] for i in 1:length(new_wt)])
                system_main.measurements[reg][:DD][m.observables[k]] =
                sum([new_wt[i] .* system_altpolicies[i].measurements[reg][:DD][m.observables[k]] for i in 1:length(new_wt)])=#
                system_main.measurements[reg][:ZZ][m.observables[k], :] =
                sum([new_wt[i] .* system_altpolicies[i].measurements[reg][:ZZ][m.observables[k], :]
                     for i in 1:(length(new_wt) - 1)]) +
                new_wt[end] .* system_taylor.measurement[:ZZ][m.observables[k], :]
                system_main.measurements[reg][:DD][m.observables[k]] =
                sum([new_wt[i] .* system_altpolicies[i].measurements[reg][:DD][m.observables[k]]
                     for i in 1:(length(new_wt) - 1)]) +
                new_wt[end] .* system_taylor.measurement[:DD][m.observables[k]]
            end
        else
            # TODO: This needs to be updated
            system_main.measurements[reg][:ZZ] .=
            sum([new_wt[i] .* system_altpolicies[i].measurements[reg][:ZZ] for i in 1:length(new_wt)])
            system_main.measurements[reg][:DD] .=
            sum([new_wt[i] .* system_altpolicies[i].measurements[reg][:DD] for i in 1:length(new_wt)])
        end

        if has_pseudo
            if has_fwd_looking_pseudo
                for k in get_setting(m, :forward_looking_pseudo_observables)
                    # COMMENTED CODE IS THE DESIRED VERSION in some form or fashion. May want to write a helper function which
                    # figures out whether a given system is regime-switching or not.
                    #=                    system_main.pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :] =
                    sum([new_wt[i] .* system_altpolicies[i].pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :]
                    for i in 1:length(new_wt)])
                    system_main.pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]] =
                    sum([new_wt[i] .* system_altpolicies[i].pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]]
                    for i in 1:length(new_wt)])=#
                    system_main.pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :] =
                    sum([new_wt[i] .* system_altpolicies[i].pseudo_measurements[reg][:ZZ_pseudo][m.pseudo_observables[k], :]
                         for i in 1:(length(new_wt) - 1)]) +
                    new_wt[end] .* system_taylor.pseudo_measurement[:ZZ_pseudo][m.pseudo_observables[k], :]
                    system_main.pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]] =
                    sum([new_wt[i] .* system_altpolicies[i].pseudo_measurements[reg][:DD_pseudo][m.pseudo_observables[k]]
                         for i in 1:(length(new_wt) - 1)]) +
                    new_wt[end] .* system_taylor.pseudo_measurement[:DD_pseudo][m.pseudo_observables[k]]
                end
            else
                # TODO: This needs to be updated
                system_main.pseudo_measurements[reg][:ZZ_pseudo] .=
                sum([new_wt[i] .* system_altpolicies[i].pseudo_measurements[reg][:ZZ_pseudo] for i in 1:length(new_wt)])
                system_main.pseudo_measurements[reg][:DD_pseudo] .=
                sum([new_wt[i] .* system_altpolicies[i].pseudo_measurements[reg][:DD_pseudo] for i in 1:length(new_wt)])
            end
        end
    end

    return system_main
end


function compute_system_helper(m::AbstractDSGEModel{T}; apply_altpolicy::Bool = false,
                               tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}

    if tvis
        @assert haskey(get_settings(m), :tvis_information_set) "The setting :tvis_information_set is not defined"
        n_tvis = haskey(get_settings(m), :tvis_regime_eqcond_info) ? length(get_setting(m, :tvis_replace_eqcond_func_dict)) : 1
        if n_tvis == 1 && haskey(get_settings(m), :tvis_regime_eqcond_info)
            if haskey(get_settings(m), :regime_eqcond_info)
                if get_setting(m, :tvis_regime_eqcond_info)[1] != get_setting(m, :replace_eqcond_func_dict)
                    warn_str = "The dictionary of functions in the Setting :regime_eqcond_info does not match the one specified " *
                    "by the length-one Setting :tvis_regime_eqcond_info. Replacing :replace_eqcond_func_dict with the dictionary " *
                    "of functions contained in :tvis_regime_eqcond_info."
                    @warn warn_str
                    m <= Setting(:regime_eqcond_info, get_setting(m, :tvis_replace_eqcond_func_dict)[1])
                end
            else
                m <= Setting(:regime_eqcond_info, get_setting(m, :tvis_replace_eqcond_func_dict)[1])
            end
        end

        if n_tvis > 1 # case of n_tvis = 1 handled below to avoid constructing redundant TimeVaryingInformationSetSystem
            @assert haskey(get_settings(m), :tvis_select_system) "The setting :tvis_select_system is not defined"
            tvis_sys = compute_tvis_system(m; apply_altpolicy = apply_altpolicy, verbose = verbose)
            transition_eqns = Transition{T}[tvis_sys[select, reg, :transition] for (reg, select) in enumerate(tvis_sys[:select])]
            return RegimeSwitchingSystem(transition_eqns, tvis_sys[:measurements], tvis_sys[:pseudo_measurements])
        end
    end

    solution_method = get_setting(m, :solution_method)

    regime_switching = haskey(get_settings(m), :regime_switching) && get_setting(m, :regime_switching)
    n_regimes        = (regime_switching && haskey(get_settings(m), :n_regimes)) ? get_setting(m, :n_regimes) : 1
    n_hist_regimes   = (regime_switching && haskey(get_settings(m), :n_hist_regimes)) ? get_setting(m, :n_hist_regimes) : 1

    # Solve model
    if regime_switching
        if solution_method == :gensys
            # Determine which regimes should use gensys2/gensys
            if haskey(get_settings(m), :gensys2) && get_setting(m, :gensys2)
                first_gensys2_regime = if haskey(get_settings(m), :regime_eqcond_info)
                    minimum(collect(keys(get_setting(m, :regime_eqcond_info))))
                else
                    GensysError("No equilibrium conditions in any regime are being temporarily replaced, " *
                                "but the setting :gensys2 is true.")
                end
                last_gensys2_regime = haskey(get_settings(m), :temporary_zlb_length) ?
                min(first_gensys2_regime + get_setting(m, :temporary_zlb_length), n_regimes) : n_regimes #NOTE removed a +1 here--if tests start failing, check here first

                gensys_regimes = UnitRange{Int}[1:(first_gensys2_regime - 1)]
                if last_gensys2_regime != n_regimes
                    append!(gensys_regimes, [(last_gensys2_regime + 1):n_regimes])
                end
                gensys2_regimes = [first_gensys2_regime-1:last_gensys2_regime]
            else
                gensys2_regimes = Vector{UnitRange{Int}}(undef, 0)
                gensys_regimes  = UnitRange{Int}[1:n_regimes]
            end

            # Solve!
            TTTs, RRRs, CCCs = solve(m; regime_switching = regime_switching,
                                     regimes = collect(1:n_regimes),
                                     gensys_regimes = gensys_regimes,
                                     gensys2_regimes = gensys2_regimes,
                                     verbose = verbose)

            transition_equations = Vector{Transition{T}}(undef, n_regimes)
            for i = 1:n_regimes
                transition_equations[i] = Transition(TTTs[i], RRRs[i], CCCs[i])
            end

            # Infer which measurement and pseudo-measurement equations to use
            type_tuple = (typeof(m), Vector{Matrix{T}}, Vector{Matrix{T}}, Vector{Vector{T}})
            #             has_pseudo = true
            has_pseudo = hasmethod(pseudo_measurement, type_tuple) ||
            hasmethod(pseudo_measurement, (typeof(m), Matrix{T}, Matrix{T}, Vector{T}))
            if tvis
                if hasmethod(measurement, type_tuple)
                    measurement_equations = measurement(m, TTTs, RRRs, CCCs;
                                                        information_set = get_setting(m, :tvis_information_set))
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
        else
            error("Regime switching with the Klein solution algorithm has not been implemented.")
        end
    else
        if solution_method == :gensys

            TTT, RRR, CCC = solve(m; verbose = verbose)
            transition_equation = Transition(TTT, RRR, CCC)

            # Solve measurement equation
            measurement_equation = measurement(m, TTT, RRR, CCC)

            #=          # TODO: delete this chunk of code, can't have an equilibrium w/ permanent ZLB, will trigger endogeneity problem
            if get_setting(m, :alternative_policy).eqcond == zero_rate_eqcond && apply_altpolicy
            measurement_equation[:DD][get_observables(m)[:obs_nominalrate]] = m[:Rstarn]
            end=#

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

"""
```
compute_system(m; apply_altpolicy = false,
            check_system = false, get_system = false,
            get_population_moments = false, use_intercept = false,
            tvis::Bool = false, verbose = :high)
            compute_system(m, data; apply_altpolicy = false,
            check_system = false, get_system = false,
            get_population_moments = false,
            tvis::Bool = false, verbose = :high)
```
Given the current model parameters, compute the DSGE-VAR or DSGE-VECM system
corresponding to model `m`. If a matrix `data` is also passed, then
the VAR is estimated on `data` using the DSGE `m` as a prior
with weight λ.

### Keyword Arguments
* `check_system::Bool`: see `?compute_system` that takes the input `m::AbstractDSGEModel`
   and `system::System`.
* `get_system::Bool`: see Outputs
* `get_population_moments::Bool`: see Outputs
* `use_intercept::Bool`: use an intercept term when computing the OLS estimate of the VAR system.
* `tvis::Bool` indicates whether the state-space system involves time-varying information sets.

### Outputs
* If `get_system = true`:
  Returns the updated `system` whose measurement matrices `ZZ`, `DD`, and `QQ` correspond
  to the VAR or VECM specified by `m`. If `m` is an `AbstractDSGEVECMModel`,
  then the `system` and the vector implied by additional cointegrating relationships
  are returned as a 2-element tuple.
* If `get_population_moments = true`:
  Returns the limit cross product matrices that describe the DSGE implied
    population moments between the observables and their lags. If `data` is
    also passed as an input, then the sample population moments are also returned.
* Otherwise:
Returns `β` and `Σ`, the coefficients and observables covariance matrix of the VAR or VECM.
If `data` is passed in, then `β` and `Σ` are estimated from the data using `m`
as a prior with weight λ. Otherwise, `β` and `Σ` comprise the VECM approximation
of the DSGE `m`.
"""
function compute_system(m::AbstractDSGEVARModel{T}; apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false, use_intercept::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T <: Real}

    regime_switching = haskey(get_settings(m), :regime_switching) ?
    get_setting(m, :regime_switching) : false
    n_regimes        = regime_switching && haskey(get_settings(m), :n_regimes) ?
    get_setting(m, :n_regimes) : 1

    dsge = get_dsge(m)
    if regime_switching
        error("Regime switching has not been implemented for a DSGEVAR yet.")
        system = compute_system(dsge; apply_altpolicy = apply_altpolicy,
                                verbose = verbose) # This `system` is really a RegimeSwitchingSystem

        systems = Vector{System{T}}(undef, n_regimes)
        for i in 1:n_regimes
            systems[i] = compute_system(dsge, System(system, i); observables = collect(keys(get_observables(m))),
                                        shocks = collect(keys(get_shocks(m))), check_system = check_system)
        end
        system = RegimeSwitchingSystem(systems) # construct a RegimeSwitchingSystem from underlying systems
    else
        system = compute_system(dsge; verbose = verbose)
        system = compute_system(dsge, system; observables = collect(keys(get_observables(m))),
                                shocks = collect(keys(get_shocks(m))), check_system = check_system)
    end

    if get_system
        return system
    elseif regime_switching
        EEs, MMs = measurement_error(m; regime_switching = regime_switching, n_regimes = n_regimes)
        out = get_population_moments ? Vector{Tuple{3, Matrix{T}}}(undef, n_regimes) :
        Vector{Tuple{2, Matrix{T}}}(undef, n_regimes)

        for i in 1:n_regimes
            out[i] = var_approx_state_space(system[i, :TTT], system[i, :RRR], system[i, :QQ],
                                            system[i, :DD], system[i, :ZZ], EEs[i], MMs[i], n_lags(m);
                                            get_population_moments = get_population_moments,
                                            use_intercept = use_intercept)
        end

        return out
    else
        EE, MM = measurement_error(m)

        return var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                      system[:DD], system[:ZZ], EE, MM, n_lags(m);
                                      get_population_moments = get_population_moments,
                                      use_intercept = use_intercept)
    end
end

function compute_system(m::AbstractDSGEVARModel{T}, data::Matrix{T};
                        apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T<:Real}

    if get_λ(m) == Inf
        # Then we just want the VAR approximation of the DSGE
        return compute_system(m; apply_altpolicy = apply_altpolicy,
                              check_system = check_system, get_system = get_system,
                              get_population_moments = get_population_moments, use_intercept = true,
                              verbose = verbose)
    else
        # Create a system using the method for DSGE-VARs and λ = ∞
        system = compute_system(m; apply_altpolicy = apply_altpolicy,
                                check_system = check_system,
                                get_system = true, use_intercept = true,
                                verbose = verbose)

        if get_system
            return system
        else
            EE, MM = measurement_error(m)

            lags = n_lags(m)
            YYYY, XXYY, XXXX =
            compute_var_population_moments(data, lags; use_intercept = true)
            out = var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                         system[:DD], system[:ZZ], EE, MM, n_lags(m);
                                         get_population_moments = true,
                                         use_intercept = true)

            if get_population_moments
                return out..., YYYY, XXYY, XXXX
            else
                # Compute prior-weighted population moments
                λ = get_λ(m)
                YYYYC = YYYY + λ .* out[1]
                XXYYC = XXYY + λ .* out[2]
                XXXXC = XXXX + λ .* out[3]

                # Draw stationary VAR system
                n_periods = size(data, 2) - lags
                β, Σ =  draw_stationary_VAR(YYYYC, XXYYC, XXXXC,
                                            convert(Int, floor(n_periods + λ * n_periods)),
                                            size(data, 1), lags)

                return β, Σ
            end
        end
    end
end

# Same functions as above but for AbstractDSGEVECMModel types
function compute_system(m::AbstractDSGEVECMModel{T}; apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false, use_intercept::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T<:Real}

    dsge = get_dsge(m)
    system = compute_system(dsge; verbose = verbose)

    # Use wrapper compute_system for AbstractDSGEVECMModel types
    # (as opposed to compute_system(m::AbstractDSGEModel, system::System; ...))
    system, DD_coint_add = compute_system(m, system; observables = collect(keys(get_observables(m))),
                                          cointegrating = collect(keys(get_cointegrating(m))),
                                          cointegrating_add = collect(keys(get_cointegrating_add(m))),
                                          shocks = collect(keys(get_shocks(m))), check_system = check_system,
                                          get_DD_coint_add = true)

    if get_system
        return system, DD_coint_add
    else
        EE, MM = measurement_error(m)

        return vecm_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                       system[:DD], system[:ZZ], EE, MM, n_observables(m),
                                       n_lags(m), n_cointegrating(m), n_cointegrating_add(m),
                                       DD_coint_add;
                                       get_population_moments = get_population_moments,
                                       use_intercept = use_intercept)
    end
end

function compute_system(m::AbstractDSGEVECMModel{T}, data::Matrix{T};
                        apply_altpolicy::Bool = false,
                        check_system::Bool = false, get_system::Bool = false,
                        get_population_moments::Bool = false,
                        tvis::Bool = false, verbose::Symbol = :high) where {T<:Real}

    if get_λ(m) == Inf
        # Then we just want the VECM approximation of the DSGE
        # with no additional cointegration
        return compute_system(m; apply_altpolicy = apply_altpolicy,
                              check_system = check_system, get_system = get_system,
                              get_population_moments = get_population_moments, use_intercept = true,
                              verbose = verbose)
    else
        # Create a system using the method for DSGE-VECMs and λ = ∞
        system, DD_coint_add = compute_system(m; apply_altpolicy = apply_altpolicy,
                                              check_system = check_system,
                                              get_system = true, use_intercept = true,
                                              verbose = verbose)

        if get_system
            return system
        else
            EE, MM = measurement_error(m)

            lags = n_lags(m)
            YYYY, XXYY, XXXX =
            compute_var_population_moments(data, lags; use_intercept = true)
            out = vecm_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                          system[:DD], system[:ZZ], EE, MM, size(data, 1),
                                          n_lags(m), n_cointegrating(m),
                                          n_cointegrating_add(m), DD_coint_add;
                                          get_population_moments = true,
                                          use_intercept = true)

            if get_population_moments
                return out..., YYYY, XXYY, XXXX
            else
                # Compute prior-weighted population moments
                λ = get_λ(m)
                YYYYC = YYYY + λ .* out[1]
                XXYYC = XXYY + λ .* out[2]
                XXXXC = XXXX + λ .* out[3]

                # Draw VECM system
                n_periods = size(data, 2) - lags
                β, Σ =  draw_VECM(YYYYC, XXYYC, XXXXC,
                                  convert(Int, n_periods + λ * n_periods),
                                  size(data, 1), lags, n_cointegrating(m))

                return β, Σ
            end
        end
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
compute_system(m::AbstractDSGEModel, system::System;
observables::Vector{Symbol} = collect(keys(m.observables)),
pseudo_observables::Vector{Symbol} = collect(keys(m.pseudo_observables)),
states::Vector{Symbol} = vcat(collect(keys(m.endogenou_states)),
collect(keys(m.endogenous_states_augmented)))
shocks::Vector{Symbol} = collect(keys(m.exogenous_shocks)),
zero_DD = false, zero_DD_pseudo = false)
compute_system(m::AbstractDSGEVECMModel, system::System;
observables::Vector{Symbol} = collect(keys(m.observables)),
pseudo_observables::Vector{Symbol} = collect(keys(m.pseudo_observables)),
cointegrating::Vector{Symbol} = collect(keys(m.cointegrating)),
states::Vector{Symbol} = vcat(collect(keys(m.endogenou_states)),
collect(keys(m.endogenous_states_augmented)))
shocks::Vector{Symbol} = collect(keys(m.exogenous_shocks)),
zero_DD = false, zero_DD_pseudo = false)
```
computes the corresponding transition and measurement
equations specified by the keywords (e.g. `states`, `pseudo_observables`)
using the existing `ZZ`, `ZZ_pseudo`, and `TTT` matrices in `system`.

Note that this function does not update the EE matrix, which is
set to all zeros. To incorporate measurement errors, the user
must specify the EE matrix after applying compute_system.

### Keywords
* `observables`: variables that should be
entered into the new `ZZ` and `DD` matrices.
The `observables` can be both Observables and PseudoObservables,
but they must be an element of the system already.
* `pseudo_observables`: variables that should be
entered into the new `ZZ_pseudo` and `DD_pseudo` matrices.
The `observables` can be both Observables and PseudoObservables,
but they must be an element of the system already.
* `cointegrating`: variables that should be
entered into the new `ZZ` and `DD` matrices as cointegrating relationships.
The `observables` can be both Observables and PseudoObservables,
but they must be an element of the system already.
* `states`: variables that should be
entered into the new `TTT` and `RRR` matrices as states.
They must be existing states.
* `shocks`: variables that should be
entered into the new `RRR` and `QQ` matrices as shocks.
They must be existing exogenous shocks.
"""
function compute_system(m::AbstractDSGEModel{S}, system::System;
                        observables::Vector{Symbol} = collect(keys(m.observables)),
                        pseudo_observables::Vector{Symbol} =
                        collect(keys(m.pseudo_observables)),
                        states::Vector{Symbol} =
                        vcat(collect(keys(m.endogenous_states)),
                             collect(keys(m.endogenous_states_augmented))),
                        shocks::Vector{Symbol} = collect(keys(m.exogenous_shocks)),
                        zero_DD::Bool = false, zero_DD_pseudo::Bool = false,
                        check_system::Bool = true) where {S<:Real}

    # Set up indices
    oid  = m.observables # observables indices dictionary
    pid  = m.pseudo_observables # pseudo observables indices dictionary
    sid  = m.endogenous_states
    said = m.endogenous_states_augmented # states augmented

    # Find shocks to keep
    shock_inds = map(k -> m.exogenous_shocks[k], shocks)
    Qout = system[:QQ][shock_inds, shock_inds]

    # Find states to keep
    if !issubset(states, vcat(collect(keys(sid)), collect(keys(said))))
        false_states = setdiff(states, vcat(collect(keys(sid)), collect(keys(said))))
        error("The following states in keyword `states` do not exist in `system`: " *
              join(string.(false_states), ", "))
    elseif !isempty(setdiff(vcat(collect(keys(sid)), collect(keys(said))), states))
        which_states = Vector{Int}(undef, length(states))
        for i = 1:length(states)
            which_states[i] = haskey(sid, states[i]) ? sid[states[i]] : said[states[i]]
        end
        Tout = system[:TTT][which_states, which_states]
        Rout = system[:RRR][which_states, shock_inds]
        Cout = system[:CCC][which_states]
    else
        which_states = 1:n_states_augmented(m)
        Tout = copy(system[:TTT])
        Rout = system[:RRR][:, shock_inds]
        Cout = copy(system[:CCC])
    end

    # Compute new ZZ and DD matrices if different observables than current system
    if !isempty(symdiff(observables, collect(keys(oid))))
        Zout = zeros(S, length(observables), size(Tout, 1))
        Dout = zeros(S, length(observables))
        for (i, obs) in enumerate(observables)
            Zout[i, :], Dout[i] = if haskey(oid, obs)
                system[:ZZ][oid[obs], which_states], zero_DD ? zero(S) : system[:DD][oid[obs]]
            elseif haskey(pid, obs)
                system[:ZZ_pseudo][pid[obs], which_states], zero_DD ? zero(S) : system[:DD_pseudo][pid[obs]]
            else
                error("Observable/PseudoObservable $obs cannot be found in the DSGE model $m")
            end
        end
    else
        Zout = copy(system[:ZZ])[:, which_states]
        Dout = zero_DD ? zeros(S, size(Zout, 1)) : copy(system[:DD])
    end

    Eout = zeros(S, length(observables), length(observables)) # measurement errors are set to zero

    # Compute new ZZ_pseudo, DD_pseudo if different pseudo_observables than current system
    if !isempty(symdiff(pseudo_observables, collect(keys(pid))))
        Zpseudoout = zeros(S, length(pseudo_observables), size(Tout, 1))
        Dpseudoout = zeros(S, length(pseudo_observables))
        for (i, pseudoobs) in enumerate(pseudo_observables)
            Zpseudoout[i, :], Dpseudoout[i] = if haskey(oid, pseudoobs)
                system[:ZZ][oid[pseudoobs], which_states], zero_DD_pseudo ?
                zero(S) : system[:DD][oid[pseudoobs]]
            elseif haskey(pid, pseudoobs)
                system[:ZZ_pseudo][pid[pseudoobs], which_states], zero_DD_pseudo ?
                zero(S) : system[:DD_pseudo][pid[pseudoobs]]
            else
                error("Observable/PseudoObservable $pseudoobs cannot be found in the DSGE model $m")
            end
        end
    else
        Zpseudoout = copy(system[:ZZ_pseudo])[:, which_states]
        Dpseudoout = zero_DD_pseudo ? zeros(S, size(Zpseudoout, 1)) : copy(system[:DD_pseudo])
    end

    if check_system
        @assert size(Zout, 2) == size(Tout, 1) "Dimension 2 of new ZZ ($(size(Zout,2))) does not match dimension of states ($(size(Tout,1)))."
        @assert size(Qout, 1) == size(Rout, 2) "Dimension 2 of new RRR ($(size(Zout,2))) does not match dimension of shocks ($(size(Qout,1)))."
    end

    # Construct new system
    return System(Transition(Tout, Rout, Cout),
                  Measurement(Zout, Dout, Qout, Eout),
                  PseudoMeasurement(Zpseudoout, Dpseudoout))
end

function compute_system(m::AbstractDSGEVECMModel{S}, system::System;
                        observables::Vector{Symbol} = collect(keys(get_observables(get_dsge(m)))),
                        cointegrating::Vector{Symbol} = Vector{Symbol}(undef, 0),
                        cointegrating_add::Vector{Symbol} = Vector{Symbol}(undef, 0),
                        pseudo_observables::Vector{Symbol} =
                        collect(keys(get_dsge(m).pseudo_observables)),
                        states::Vector{Symbol} =
                        vcat(collect(keys(get_dsge(m).endogenous_states)),
                             collect(keys(get_dsge(m).endogenous_states_augmented))),
                        shocks::Vector{Symbol} = collect(keys(get_dsge(m).exogenous_shocks)),
                        zero_DD::Bool = false, zero_DD_pseudo::Bool = false,
                        get_DD_coint_add::Bool = false,
                        check_system::Bool = true) where {S<:Real}
    # Cointegrating relationships should exist as observables/pseudo_observables already
    # in the underlying DSGE. We assume cointegrating relationships come after normal observables.
    # Default behavior is to recreate the underlying DSGE's state space representation, however.
    sys = compute_system(get_dsge(m), system; observables = vcat(observables, cointegrating),
                         pseudo_observables = pseudo_observables,
                         states = states, shocks = shocks, zero_DD = zero_DD,
                         zero_DD_pseudo = zero_DD_pseudo, check_system = check_system)
    if get_DD_coint_add
        mtype = typeof(m)
        DD_coint_add = if hasmethod(compute_DD_coint_add, (mtype, Vector{Symbol}))
            compute_DD_coint_add(m, cointegrating_add)
        elseif hasmethod(compute_DD_coint_add, (mtype, ))
            compute_DD_coint_add(m)
        else
            compute_DD_coint_add(m, sys, cointegrating_add)
        end
        return sys, DD_coint_add
    else
        return sys
    end
end

"""
```
function compute_DD_coint_add(m::AbstractDSGEVECMModel{S}, system::System,
cointegrating_add::Vector{Symbol}) where {S <: Real}
```
computes `DD_coint_add` for a `DSGEVECM` model. This vector
holds additional cointegrating relationships that do not require
changes to the `ZZ` matrix.

### Note
We recommend overloading this function if there are
cointegrating relationships which a user does not want
to add to the underlying DSGE. The function `compute_system`
first checks for a method `compute_DD_coint_add` that takes
a Type tuple of `(model_type, Vector{Symbol})` and then `(model_type, )`
before calling this method.

This function is generally intended to be internal. As an example of
other such functions, `eqcond` must be user-defined but
is primarily used internally and not directly called by the user in a script.
"""
function compute_DD_coint_add(m::AbstractDSGEVECMModel{S}, system::System,
                              cointegrating_add::Vector{Symbol}) where {S <: Real}
    if !isempty(cointegrating_add)
        oid = get_observables(get_dsge(m))
        pid = get_pseudo_observables(get_dsge(m))
        Dout = zeros(S, length(cointegrating_add))
        for (i, obs) in enumerate(cointegrating_add)
            Dout[i] = if haskey(oid, obs)
                system[:DD][oid[obs]]
            elseif haskey(pid, obs)
                system[:DD_pseudo][pid[obs]]
            else
                error("Observable/PseudoObservable $obs cannot be found in the DSGE model $m")
            end
        end
        return Dout
    else
        @warn "No additional cointegrating relationships specified. Returning an empty vector."
        return Vector{S}(undef, 0)
    end
end

"""
```
compute_system_function(system::System{S}) where S<:AbstractFloat
```

### Inputs

- `system::System`: The output of compute_system(m), i.e. the matrix outputs from solving a given model, m.

### Outputs

- `Φ::Function`: transition equation
- `Ψ::Function`: measurement equation
- `F_ϵ::Distributions.MvNormal`: shock distribution
- `F_u::Distributions.MvNormal`: measurement error distribution
"""
function compute_system_function(system::System{S}) where S<:AbstractFloat
    # Unpack system
    TTT    = system[:TTT]
    RRR    = system[:RRR]
    CCC    = system[:CCC]
    QQ     = system[:QQ]
    ZZ     = system[:ZZ]
    DD     = system[:DD]
    EE     = system[:EE]

    # Define transition and measurement functions
    @inline Φ(s_t1::Vector{S}, ϵ_t::Vector{S}) = TTT*s_t1 + RRR*ϵ_t + CCC
    @inline Ψ(s_t::Vector{S}) = ZZ*s_t + DD

    # Define shock and measurement error distributions
    nshocks = size(QQ, 1)
    nobs    = size(EE, 1)
    F_ϵ = Distributions.MvNormal(zeros(nshocks), QQ)
    F_u = Distributions.MvNormal(zeros(nobs),    EE)

    return Φ, Ψ, F_ϵ, F_u
end

function zero_system_constants(system::System{S}) where S<:AbstractFloat
    system = copy(system)

    system.transition.CCC = zeros(size(system[:CCC]))
    system.measurement.DD = zeros(size(system[:DD]))
    system.pseudo_measurement.DD_pseudo = zeros(size(system[:DD_pseudo]))

    return system
end

function zero_system_constants(system::RegimeSwitchingSystem{S}) where S<:AbstractFloat
    system = copy(system)

    for i in 1:n_regimes(system)
        system.transitions[i].CCC = zeros(size(system[i, :CCC]))
        system.measurements[i].DD = zeros(size(system[i, :DD]))
        system.pseudo_measurements[i].DD_pseudo = zeros(size(system[i, :DD_pseudo]))
    end

    return system
end

"""
```
var_approx_state_space(TTT, RRR, QQQ, DD, ZZ, EE, MM, p; get_population_moments = false,
use_intercept = false) where {S<:Real}
```
computes the VAR(p) approximation of the linear state space system

```
sₜ = TTT * sₜ₋₁ + RRR * ϵₜ,
yₜ = ZZ * sₜ + DD + uₜ,
```
where the disturbances are assumed to follow
```
ϵₜ ∼ 𝒩 (0, QQ),
uₜ = ηₜ + MM * ϵₜ,
ηₜ ∼ 𝒩 (0, EE).
```
The `MM` matrix implies
```
cov(ϵₜ, uₜ) = QQ * MM'.
```

### Outputs
If `get_population_moments = false`:
* `β`: VAR(p) coefficients
* `Σ`: innovations variance-covariance matrix for the VAR(p) representation
```
yₜ = Xₜβ + μₜ
```
where `Xₜ` appropriately stacks the constants and `p` lags of `yₜ`, and `μₜ ∼ 𝒩 (0, Σ)`.

If `get_population_moments = true`: we return the limit cross product matrices.
* `yyyyd`: 𝔼[y,y]
* `XXXXd`: 𝔼[y,X(lag rr)]
* `XXyyd`: 𝔼[X(lag rr),X(lag ll)]

Using these matrices, the VAR(p) representation is given by
```
β = XXXXd \\ XXyyd
Σ = yyyyd - XXyyd' * β
```

The keyword `use_intercept` specifies whether or not to use an
intercept term in the VAR approximation.
"""
function var_approx_state_space(TTT::AbstractMatrix{S}, RRR::AbstractMatrix{S},
                                QQ::AbstractMatrix{S}, DD::AbstractVector{S},
                                ZZ::AbstractMatrix{S}, EE::AbstractMatrix{S},
                                MM::AbstractMatrix{S}, p::Int;
                                get_population_moments::Bool = false,
                                use_intercept::Bool = false) where {S<:Real}

    n_obs = size(ZZ, 1)

    HH = EE + MM * QQ * MM'
    VV = QQ * MM'

    ## Compute p autocovariances

    ## Initialize Autocovariances
    GAMM0 = zeros(S, n_obs ^ 2, p + 1)
    GA0 =  solve_discrete_lyapunov(TTT, RRR * QQ * RRR')
    Gl   = ZZ * GA0 * ZZ' + ZZ * RRR * VV + (ZZ * RRR * VV)' + HH
    GAMM0[:, 1] = vec(Gl)
    TTl = copy(TTT)
    GA0ZZ = GA0 * ZZ'
    RRRVV = RRR * VV
    for l = 1:p
        Gl = ZZ * TTl * (GA0ZZ + RRRVV) # ZZ * (TTl * GA0Z) * ZZ' + ZZ * (TTl * RRR * VV)
        GAMM0[:, l+1] = vec(Gl)
        TTl = TTl * TTT
    end

    ## Create limit cross product matrices
    yyyyd = zeros(S, n_obs, n_obs)
    if use_intercept
        XXXXd = zeros(S, 1 + p * n_obs, 1 + p * n_obs)
        yyXXd = zeros(S, n_obs, 1 + p * n_obs)

        XXXXd[1, 1] = 1.
        XXXXd[1, 2:1 + p * n_obs] = repeat(DD', 1, p) # same as kron(ones(1, p), DD')
        XXXXd[2:1 + p * n_obs, 1] = repeat(DD, p)     # same as kron(ones(p), DD)
        yyXXd[:, 1] = DD
    else
        XXXXd = zeros(S, p * n_obs, p * n_obs)
        yyXXd = zeros(S, n_obs, p * n_obs)
    end
    DDDD = DD * DD'
    yyyyd = reshape(GAMM0[:, 1], n_obs, n_obs) + DDDD

    shift = use_intercept ? 1 : 0 # for constructing XXXXd, to select the right indices
    for rr = 1:p
        # 𝔼[yy,x(lag rr)]
        yyXXd[:, n_obs * (rr - 1) + 1 + shift:n_obs * rr + shift] =
        reshape(GAMM0[:, rr + 1], n_obs, n_obs) + DDDD

        # 𝔼[x(lag rr),x(lag ll)]
        for ll = rr:p
        yyyydrrll = reshape(GAMM0[:, ll - rr + 1], n_obs, n_obs) + DDDD
        XXXXd[n_obs * (rr - 1) + 1 + shift:n_obs * rr + shift,
              n_obs * (ll - 1) + 1 + shift:n_obs * ll + shift] = yyyydrrll
        XXXXd[n_obs * (ll - 1) + 1 + shift:n_obs * ll + shift,
              n_obs * (rr - 1) + 1 + shift:n_obs * rr + shift] = yyyydrrll'
        end
    end

    XXyyd = convert(Matrix{S}, yyXXd')

    if get_population_moments
        return yyyyd, XXyyd, XXXXd
    else
        β = \(XXXXd, XXyyd)
        Σ = yyyyd - XXyyd' * β
        return β, Σ
    end
end

"""
```
vecm_approx_state_space(TTT, RRR, QQQ, DD, ZZ, EE, MM, n_obs, p, n_coint,
n_coint, n_coint_add, DD_coint_add; get_population_moments = false,
use_intercept = false) where {S<:Real}
```
computes the VECM(p) approximation of the linear state space system

```
sₜ = TTT * sₜ₋₁ + RRR * ϵₜ,
yₜ = ZZ * sₜ + DD + uₜ,
```
where the disturbances are assumed to follow
```
ϵₜ ∼ 𝒩 (0, QQ),
uₜ = ηₜ + MM * ϵₜ,
ηₜ ∼ 𝒩 (0, EE).
```
The `MM` matrix implies
```
cov(ϵₜ, uₜ) = QQ * MM'.
```

### Outputs
If `get_population_moments = false`:
* `β`: VECM(p) coefficients. The first `n_coint + n_coint_add`
coefficients for each observable comprise the error correction terms,
while the following `1 + p * n_obs` terms are the VAR terms.
* `Σ`: innovations variance-covariance matrix for the VECM(p) representation
```
Δyₜ = eₜβₑ + Xₜβᵥ + μₜ
```
where `βₑ` are the coefficients for the error correction terms;
`eₜ` are the error correction terms specifying the cointegrating relationships;
`βᵥ` are the coefficients for the VAR terms; `Xₜ` appropriately stacks
the constants and `p` lags of `Δyₜ`; and `μₜ ∼ 𝒩 (0, Σ)`.

Note that the error correction terms satisfy the mapping
`eₜ' = C * yₜ₋₁`, where `C` is a matrix.

If `get_population_moments = true`: we return the limit cross product matrices.
* `yyyyd`: 𝔼[y,y]
* `XXXXd`: 𝔼[y,X(lag rr)]
* `XXyyd`: 𝔼[X(lag rr),X(lag ll)]

Note that in the rows of `XXyyd` and the rows and columns of `XXXXd`,
the cointegrating relationships are stacked above the constants and
lagged `Δyₜ`.

Using these matrices, the VAR(p) representation is given by
```
β = XXXXd \\ XXyyd
Σ = yyyyd - XXyyd' * β,
```
where `β` has dimensions `n_obs × (n_coint + n_coint_add + 1 + p * n_obs)`,
and `Σ` has dimensions `n_obs × n_obs`.

The keyword `use_intercept` specifies whether or not to use an
intercept term in the VECM approximation.
"""
function vecm_approx_state_space(TTT::AbstractMatrix{S}, RRR::AbstractMatrix{S},
                                 QQ::AbstractMatrix{S}, DD::AbstractVector{S},
                                 ZZ::AbstractMatrix{S}, EE::AbstractMatrix{S},
                                 MM::AbstractMatrix{S}, n_obs::Int, p::Int,
                                 n_coint::Int, n_coint_add::Int = 0,
                                 DD_coint_add::AbstractVector{S} = Vector{S}(undef, 0);
                                 get_population_moments::Bool = false,
                                 use_intercept::Bool = false,
                                 test_GA0::AbstractMatrix{S} =
                                 Matrix{S}(undef, 0, 0)) where {S <: Real}

    # n_obs is number of observables, n_coint is number of cointegrating relationships
    n_coint_all = n_coint + n_coint_add

    # Create variance-covariance matrices w/measurement error included
    HH = EE + MM * QQ * MM'
    VV = QQ * MM'

    ## Compute p autocovariances

    ## Initialize Autocovariances
    GAMM0 = zeros(S, (n_obs + n_coint) ^ 2, p + 1)
    GA0 = isempty(test_GA0) ? solve_discrete_lyapunov(TTT, RRR * QQ * RRR') : test_GA0 # Matlab code uses a different numerical procedure -> some difference in this matrix
    Gl   = ZZ * GA0 * ZZ' + ZZ * RRR * VV + (ZZ * RRR * VV)' + HH
    GAMM0[:, 1] = vec(Gl)
    TTl = copy(TTT)
    GA0ZZ = GA0 * ZZ'
    RRRVV = RRR * VV
    for l = 1:p
        Gl = ZZ * TTl * (GA0ZZ + RRRVV) # ZZ * (TTl * GA0) * ZZ' + ZZ * (TTl * RRR * VV)
        GAMM0[:, l+1] = vec(Gl)
        TTl = TTl * TTT
    end

    ## Create limit cross product matrices
    DDDD = DD * DD'
    yyyyd_coint0 = reshape(GAMM0[:, 1], n_obs + n_coint, n_obs + n_coint) + DDDD
    yyyyd_coint1 = reshape(GAMM0[:, 2], n_obs + n_coint, n_obs + n_coint) + DDDD
    yyyyd = yyyyd_coint0[1:n_obs, 1:n_obs]

    # n_coint_add are treated as the first set of variables in XX
    # n_coint    are treated as the second set of variables in XX
    # composition: n_coint_add - n_coint - constant - lags
    if use_intercept
        XXXXd = zeros(S, 1 + p * n_obs + n_coint_all, 1 + p * n_obs + n_coint_all)
        yyXXd = zeros(S, n_obs, 1 + p * n_obs + n_coint_all)

        # 𝔼[x(n_coint), x(n_coint)]
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_add + 1:n_coint_all] = yyyyd_coint0[n_obs + 1:n_obs + n_coint, n_obs + 1:n_obs + n_coint]

        # 𝔼[x(const), x(const)]
        XXXXd[n_coint_all + 1, n_coint_all + 1] = 1.

        # 𝔼[x(n_coint), x(const)]
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_all + 1] = DD[n_obs + 1:n_obs + n_coint]
        XXXXd[n_coint_all + 1, n_coint_add + 1:n_coint_all] = DD[n_obs + 1:n_obs + n_coint]

        # 𝔼[x(const), x(lags)]
        XXXXd[n_coint_all + 1, n_coint_all + 2:n_coint_all + 1 + p * n_obs] = repeat(DD[1:n_obs]', 1, p)
        XXXXd[n_coint_all + 2:n_coint_all + 1 + p * n_obs, n_coint_all + 1] =
        XXXXd[n_coint_all + 1, n_coint_all + 2:n_coint_all + 1 + p * n_obs]' # same as kron(ones(p), DD[1:n_obs]) but avoids calculation

        # 𝔼[yy, x(n_coint)]
        yyXXd[:, n_coint_add + 1:n_coint_all] = yyyyd_coint1[1:n_obs, n_obs + 1:n_obs + n_coint]

        # 𝔼[yy, x(const)]
        yyXXd[:, n_coint_all + 1] = DD[1:n_obs]
    else
        XXXXd = zeros(S, p * n_obs + n_coint_all, p * n_obs + n_coint_all)
        yyXXd = zeros(S, n_obs, p * n_obs + n_coint_all)

        # 𝔼[x(n_coint), x(n_coint)]
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_add + 1:n_coint_all] =
        yyyyd_coint0[n_obs + 1:n_obs + n_coint, n_obs + 1:n_obs + n_coint]

        # 𝔼[yy, x(n_coint)]
        yyXXd[:, n_coint_add + 1:n_coint_all] = yyyyd_coint1[1:n_obs, n_obs + 1:n_obs + n_coint]
    end

    if n_coint_add > 0
        DD_coint_add_div2 = DD_coint_add ./ 2
        if use_intercept
            # 𝔼[yy, x(n_coint_add)]
            yyXXd[:, 1:n_coint_add] = DD[1:n_obs] * DD_coint_add_div2

            # 𝔼[x(n_coint_add), x(n_coint_add)]
            XXXXd[1:n_coint_add, 1:n_coint_add] = DD_coint_add * DD_coint_add' ./ 3

            # 𝔼[x(n_coint_add), x(n_coint)]
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all] =
            DD_coint_add_div2 * DD[n_obs + 1: n_obs + n_coint]'
            XXXXd[1 + n_coint_add:n_coint_all, 1:n_coint_add] =
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all]' # transpose of the previous line

            # 𝔼[x(n_coint_add), x(const)]
            XXXXd[1:n_coint_add, n_coint_all + 1] = DD_coint_add_div2
            XXXXd[n_coint_all + 1, 1:n_coint_add] = DD_coint_add_div2'
        else
            # 𝔼[yy, x(n_coint_add)]
            yyXXd[:, 1:n_coint_add] = DD[1:n_obs] * DD_coint_add_div2

            # 𝔼[x(n_coint_add), x(n_coint_add)]
            XXXXd[1:n_coint_add, 1:n_coint_add] = DD_coint_add * DD_coint_add' ./ 3

            # 𝔼[x(n_coint_add), x(n_coint)]
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all] =
            DD_coint_add_div2 * DD[n_obs + 1: n_obs + n_coint]'
            XXXXd[1 + n_coint_add:n_coint_all, 1:n_coint_add] =
            XXXXd[1:n_coint_add, 1 + n_coint_add:n_coint_all]' # transpose of the previous line
        end
    end

    shift = use_intercept ? 1 : 0 # for constructing XXXXd, to select the right indices
    for rr = 1:p
        # 𝔼[yy, x(lag rr)]
        yyyyd_cointrr = reshape(GAMM0[:, rr + 1], n_obs + n_coint, n_obs + n_coint) + DDDD
        yyXXd[:, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] =
        yyyyd_cointrr[1:n_obs, 1:n_obs]

        if n_coint_add > 0
            # 𝔼[x(n_coint_add), x(lag rr)]
            XXXXd[1:n_coint_add, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] =
            DD_coint_add_div2 * DD[1:n_obs]'
            XXXXd[n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift, 1:n_coint_add] =
            XXXXd[1:n_coint_add, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift]'
        end

        # 𝔼[x(n_coint), x(lag rr)]
        yyyyd_cointrr1 = reshape(GAMM0[:, rr], n_obs + n_coint, n_obs + n_coint) + DDDD
        XXXXd[n_coint_add + 1:n_coint_all, n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] =
        yyyyd_cointrr1[n_obs + 1:n_obs + n_coint, 1:n_obs]
        XXXXd[n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift, n_coint_add + 1:n_coint_all] =
        yyyyd_cointrr1[n_obs + 1:n_obs + n_coint, 1:n_obs]'

        # 𝔼[x(lag rr), x(lag ll)]
        for ll = rr:p
            yyyyd_cointrrll = reshape(GAMM0[:, ll - rr + 1], n_obs + n_coint, n_obs + n_coint) + DDDD
            XXXXd[n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift,
                  n_coint_all + 1 + n_obs * (ll - 1) + shift:n_coint_all + n_obs * ll + shift] = yyyyd_cointrrll[1:n_obs, 1:n_obs]
            XXXXd[n_coint_all + 1 + n_obs * (ll - 1) + shift:n_coint_all + n_obs * ll + shift,
                  n_coint_all + 1 + n_obs * (rr - 1) + shift:n_coint_all + n_obs * rr + shift] = yyyyd_cointrrll[1:n_obs, 1:n_obs]'
        end
    end

    XXyyd = convert(Matrix{S}, yyXXd')

    if get_population_moments
        return yyyyd, XXyyd, XXXXd
    else
        β = XXXXd \ XXyyd
        Σ = yyyyd - XXyyd' * β
        return β, Σ
    end
end

"""
```
k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, k, permanent_t = length(TTTs))
```

calculates the matrices associated with the expected state `k` periods ahead from `t`.
This function should NOT be used with linear state space system matrices with any unit roots.

The `TTT` and `CCC` inputs are the transition matrix and constant vector associated with
the current period `t`, while the `TTTs` and `CCCs` are vectors containing the time-varying
transition matrices and constant vectors, such that `TTTs[t]` retrieves the time-varying
transition matrix associated with period `t` and `TTTs[t + k]` retrieves the time-varying
transition matrix associated with period `t + k`. The optional argument `permanent_t`
indicates the period for which the matrices/vectors are no longer time-varying, i.e.
if `t >= permanent_t`, then `TTTs[permanent_t]` is the transition matrix.

The formula implemented by this function is
```
𝔼ₜ[sₜ₊ₖ] = (∏ⱼ=₁ᵏ Tₜ₊ⱼ) sₜ + (∑ₘ₌₁ᵏ⁻¹ (∏ⱼ₌ₘ₊₁ᵏ Tₜ₊ⱼ) Cₜ₊ₘ) + Cₜ₊ₖ.
```
Additional simplifications are made if it is known that `t + k > permanent_t`
since this implies some matrices are the same. This recognition reduces
unnecessary computations.
"""
function k_periods_ahead_expectations(TTT::AbstractMatrix, CCC::AbstractVector,
                                      TTTs::Vector{<: AbstractMatrix}, CCCs::Vector{<: AbstractVector},
                                      t::Int, k::Int, permanent_t::Int = length(TTTs);
                                      integ_series::Bool = false)

    if isempty(TTTs) || isempty(CCCs)
        Tᵏ = TTT^k
        if all(CCC .≈ 0.)
            return Tᵏ, CCC
        else
            if integ_series
                T_memo = Dict{Int, typeof(TTT)}()
                T_memo[1] = TTT
                for i in 2:k
                    T_memo[i] = T_memo[i - 1] * TTT
                end
                Tᵏsum = sum([T_memo[j] for j in 1:k])
            else
                Tᵏsum = (I - TTT) \ (I - Tᵏ)
            end
            return Tᵏ, Tᵏsum * CCC
        end
    else
        if t + k <= permanent_t
            # Cannot save computation speed by not calculating further times b/c always time-varying
            T_memo = Dict{Int, eltype(TTTs)}()
            T_memo[k] = TTTs[t + k]
            for i in (k-1):-1:1
                T_memo[i] = T_memo[i + 1] * TTTs[t + i]
            end

            C_accum = deepcopy(CCCs[t + k])
            for i in 1:(k - 1)
                C_accum .+= T_memo[i + 1] * CCCs[t + i]
            end

            return T_memo[1], C_accum
        elseif t == permanent_t
            # None of the matrices are time-varying anymore
            if all(CCCs[permanent_t] .≈ 0.)
                return TTTs[permanent_t]^k, CCCs[permanent_t]
            else
                Tᵏₜ₊₁ = TTTs[permanent_t]^k
                if integ_series
                    T_memo = Dict{Int, eltype(TTTs)}()
                    T_memo[1] = TTTs[permanent_t]
                    for i in 2:k
                        T_memo[i] = T_memo[i - 1] * TTTs[permanent_t]
                    end
                    Tᵏsum = sum([T_memo[j] for j in 1:k])
                else
                    Tᵏsum = (I - TTTs[permanent_t]) \ (I - Tᵏₜ₊₁)
                end
                return Tᵏₜ₊₁, Tᵏsum * CCCs[permanent_t]
            end
        else
            # Computation time can be saved by realizing some matrices are not time-varying
            h = (permanent_t - 1) - t # last time of time-variation is permanent_t - 1
            Tᵏ⁻ʰₜ₊ₕ₊₁ = (k == h) ? Diagonal(ones(length(CCCs[permanent_t]))) : TTTs[permanent_t]^(k - h)

            T_memo = Dict{Int, eltype(TTTs)}()
            if h > 0
                T_memo[h] = TTTs[t + h] # maps i to ∏ⱼ₌ᵢʰ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ, T_memo[h-1] = Tₜ₊ₕ * Tₜ₊ₕ₋₁, ...
                for i in (h-1):-1:1
                    T_memo[i] = T_memo[i + 1] * TTTs[t + i]
                end
                T_accum = Tᵏ⁻ʰₜ₊ₕ₊₁ * T_memo[1]
            else
                T_accum = Tᵏ⁻ʰₜ₊ₕ₊₁
            end

            if h == 0 # Nothing to accumulate from the past
                C_accum = zeros(size(T_accum, 1))
            else
                C_accum = deepcopy(CCCs[t + h])
                for i in 1:(h - 1)
                    C_accum .+= T_memo[i + 1] * CCCs[t + i]
                end
                C_accum .= Tᵏ⁻ʰₜ₊ₕ₊₁ * C_accum
            end

            if all(CCCs[permanent_t] .≈ 0.)
                return T_accum, C_accum
            else
                if integ_series
                    T_memo = Dict{Int, eltype(TTTs)}()
                    T_memo[1] = TTTs[permanent_t]
                    for i in 2:(k - h - 1)
                        T_memo[i] = T_memo[i - 1] * TTTs[permanent_t]
                    end
                    Tᵏ⁻ᵐsum = sum([T_memo[j] for j in 1:(k - h - 1)]) + I
                else
                    Tᵏ⁻ᵐsum = (I - TTTs[permanent_t]) \ (I - Tᵏ⁻ʰₜ₊ₕ₊₁)
                end
                C_accum .+= Tᵏ⁻ᵐsum * CCCs[permanent_t]
                return T_accum, C_accum
            end
        end
    end
end

"""
```
k_periods_ahead_expected_sums(TTT, CCC, TTTs, CCCs, t, k, permanent_t = length(TTTs))
```

calculates the matrices associated with the sum of the expected states between periods
`t + 1` and  `t + k`. This function should NOT be used with
linear state space system matrices with any unit roots.

The `TTT` and `CCC` inputs are the transition matrix and constant vector associated with
the current period `t`, while the `TTTs` and `CCCs` are vectors containing the time-varying
transition matrices and constant vectors, such that `TTTs[t]` retrieves the time-varying
transition matrix associated with period `t` and `TTTs[t + k]` retrieves the time-varying
transition matrix associated with period `t + k`. The optional argument `permanent_t`
indicates the period for which the matrices/vectors are no longer time-varying, i.e.
if `t >= permanent_t`, then `TTTs[permanent_t]` is the transition matrix.

The formula implemented by this function is
```
∑ⱼ₌₁ᵏ 𝔼ₜ[sₜ₊ⱼ] = ∑ⱼ₌₁ᵏ(∏ⱼ=₁ᵏ Tₜ₊ⱼ) sₜ + ∑ᵣ₌₁ᵏ⁻¹(I + ∑ⱼ₌ᵣ₊₁ᵏ (∏ₘ₌ᵣ₊₁ʲ Tₜ₊ₘ))Cₜ₊ᵣ + Cₜ₊ₖ.
```
Additional simplifications are made if it is known that `t + k > permanent_t`
since this implies some matrices are the same. This recognition reduces
unnecessary computations.
"""
function k_periods_ahead_expected_sums(TTT::AbstractMatrix, CCC::AbstractVector,
                                       TTTs::Vector{<: AbstractMatrix}, CCCs::Vector{<: AbstractVector},
                                       t::Int, k::Int, permanent_t::Int = length(TTTs);
                                       integ_series::Bool = false)

    if integ_series # Do this by directly summing the k-periods ahead expectations. Not fully efficient but also not the typical case
        T_accum = Vector{eltype(TTTs)}(undef, k)
        C_accum = Vector{eltype(CCCs)}(undef, k)
        for i in 1:k
            T_accum[i], C_accum[i] = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i,
                                                                  integ_series = integ_series)
        end
        return sum(T_accum), sum(C_accum)
    end

    if isempty(TTTs) || isempty(CCCs)
        Tᵏsum = (I - TTT) \ (TTT - TTT^(k + 1))
        if all(CCC .≈ 0.)
            return Tᵏsum, CCC
        else
            TTTʲmemo = Vector{typeof(TTT)}(undef, k)
            TTTʲmemo[1] = TTT
            for q in 2:k
                TTTʲmemo[q] = TTTʲmemo[q - 1] * TTT
            end
            return Tᵏsum, (I + sum([(I - TTT) \ (I - TTTʲmemo[k - q + 1]) for q in 1:(k - 1)])) * CCC
        end
    else
        if t + k <= permanent_t
            # Cannot save computation speed by not calculating further times b/c always time-varying
            total_Tsum = zeros(eltype(TTTs[t]), size(TTTs[t]))
            total_Csum = zeros(eltype(CCCs[t]), size(CCCs[t]))

            for i in 1:k
                tmp1, tmp2 = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i, permanent_t,
                                                          integ_series = integ_series)
                total_Tsum .+= tmp1
                total_Csum .+= tmp2
            end

            return total_Tsum, total_Csum
        elseif t == permanent_t
            # None of the matrices are time-varying anymore
            Tᵏsum = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t]^(k + 1))
            if all(CCC .≈ 0.)
                return Tᵏsum, CCCs[permanent_t]
            else
                TTTʲmemo = Vector{eltype(TTTs)}(undef, k)
                TTTʲmemo[1] = TTTs[permanent_t]
                for q in 2:k
                    TTTʲmemo[q] = TTTʲmemo[q - 1] * TTTs[permanent_t]
                end
                return Tᵏsum, (I + sum([(I - TTTs[permanent_t]) \ (I - TTTʲmemo[k - q + 1]) for q in 1:(k - 1)])) * CCCs[permanent_t]
            end
        else
            # Computation time can be saved by realizing some matrices are not time-varying
            T_accum = Vector{eltype(TTTs)}(undef, k)
            C_accum = Vector{eltype(CCCs)}(undef, k)
            for i in 1:k
            T_accum[i], C_accum[i] = k_periods_ahead_expectations(TTT, CCC, TTTs, CCCs, t, i,
                                                                  integ_series = integ_series)
        end
        return sum(T_accum), sum(C_accum)

        # This code seems to work (match directly adding each k_period_ahead_expectations)
        # but sometimes causes mild numerical differences,
        # so we calculate the expected sum by adding each k_period_ahead_expectations.
        #=            h = (permanent_t - 1) - t # last time of time-variation is permanent_t - 1

        Tₜ₊ₕ₊₁_memo = Vector{eltype(TTTs)}(undef, k - h)
        Tₜ₊ₕ₊₁_memo[1] = TTTs[permanent_t] # maps j to Tₜ₊ₕ₊₁ʲ⁻ʰ for j in (h + 1):k or Tₜ₊ₕ₊₁ʲ for j in 1:(k-h)
        for j in 2:(k - h)
        Tₜ₊ₕ₊₁_memo[j] = TTTs[permanent_t] * Tₜ₊ₕ₊₁_memo[j - 1]
        end

        Tₜ₊ₕ₊₁ʲsum = (I - TTTs[permanent_t]) \ (TTTs[permanent_t] - TTTs[permanent_t] ^ (k - h + 1)) # ∑ⱼ₌₁ᵏ⁻ʰ Tₜ₊ₕ₊₁ʲ

        TC_memo = Vector{eltype(TTTs)}(undef, h) # matrices to be used for both accumulated TTT and CCC
        TC_memo[h] = TTTs[t + h] # maps i to ∏ⱼ₌ᵢʰ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ, T_memo[h-1] = Tₜ₊ₕ * Tₜ₊ₕ₋₁, ...
        for i in (h-1):-1:1
        TC_memo[i] = TC_memo[i + 1] * TTTs[t + i]
        end

        T1_term = Vector{eltype(TTTs)}(undef, h) # first part of the accumulated TTT matrix
        T1_term[1] = TTTs[t + 1] # maps i to ∏ⱼ₌₁ⁱ Tₜ₊ⱼ, so T_memo[h] = Tₜ₊ₕ * ⋯  * Tₜ₊₁
        for i in 2:(h - 1)
        T1_term[i] = TTTs[t + i] * T1_term[i - 1]
        end
        T1_term[h] = TC_memo[1] # This one was calculated already

        # second part of the accumulated TTT matrix
        # maps j to Tₜ₊ₕ₊₁ʲ ∏ᵣ₌₁ʰ Tₜ₊ᵣ, so T_memo[h] = Tₜ₊ₕ * ⋯  * Tₜ₊₁
        T2_term = [Tₜ₊ₕ₊₁_memo[j] * T1_term[h] for j in 1:(k - h)]

        T_accum = sum(T1_term) + sum(T2_term) # calculated accumulated matrix

        # Calculate final portions of accumulated CCC vector first
        I₊Tʲ = (I + Tₜ₊ₕ₊₁ʲsum)
        C_accum   = I₊Tʲ * CCCs[t + h]
        if any(.!(CCCs[permanent_t] .≈ 0.))
        C_accum .+= (I + (k - h - 1) .* I₊Tʲ) * CCCs[permanent_t]
        end
        if h > 1
        C1_term = Vector{eltype(TTTs)}(undef, h - 1)
        C2_term = Vector{eltype(TTTs)}(undef, h - 1)

        for q in 1:(h - 1)
        C1_term[q] = sum([prod([TTTs[t + m] for m in (q + 1):j]) for j in (q + 1):h])
        C2_term[q] = sum([Tₜ₊ₕ₊₁_memo[j] * TC_memo[q + 1] for j in 1:(k - h)])
        end
        C_accum .+= sum([(I + C1_term[q] + C2_term[q]) * CCCs[t + q] for q in 1:(h - 1)])
        end=#

        return T_accum, C_accum
    end
end
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
- `:tvis_regime_eqcond_info`: Vector of `replace_eqcond_func_dict`
to specify different sets of equilibrium conditions, which
generate different state space systems (usually).
- `:tvis_select_system`: Vector of `Int` to specify which
state space system to use when calculating the measurement
and pseudo measurement equations. These state space systems
correspond to different sets of equilibrium conditions (usually).
"""
function compute_tvis_system(m::AbstractDSGEModel{T}; apply_altpolicy::Bool = false,
                         verbose::Symbol = :high) where {T <: Real}

    # TODO: update this compute_tvis_system to compute the average over forward-looking variables
    #       when using uncertain altpol and uncertain zlb. This means that, for each
    #       vector of transition equations, we need to compute the associated average of measurement equations,
    #       which requires us to calculate the perfect credibility solutions for each vector of transition equations
    #       that are meant to have uncertainty in them.

    @assert get_setting(m, :solution_method) ==
    :gensys "Currently, the solution method must be :gensys to calculate a state-space system with time-varying information sets"

    @assert get_setting(m, :replace_eqcond) "The setting :replace_eqcond must be true to calculate a state-space system with time-varying information sets"

    # :regime_dates should have the same number of possible regimes. Any differences in eqcond
    # should be specified by tvis_regime_eqcond_info
    tvis_infoset                  = get_setting(m, :tvis_information_set)
    tvis_regime_eqcond_info = get_setting(m, :tvis_replace_eqcond_func_dict)
    tvis_select                   = get_setting(m, :tvis_select_system)
    regime_switching              = get_setting(m, :regime_switching)

    n_tvis           = length(tvis_replace_eqcond_func_dict)
    n_regimes        = regime_switching && haskey(get_settings(m), :n_regimes) ?
    get_setting(m, :n_regimes) : 1
    n_hist_regimes   = regime_switching && haskey(get_settings(m), :n_hist_regimes) ?
    get_setting(m, :n_hist_regimes) : 1

    # determine which regimes use gensys/gensys2
    first_gensys2_regime = haskey(get_settings(m), :replace_eqcond_func_dict) ? min(collect(keys(get_setting(m, :replace_eqcond_func_dict)))...) : n_hist_regimes + 1
    last_gensys2_regime = haskey(get_settings(m), :temporary_zlb_length) ? first_gensys2_regime + get_setting(m, :temporary_zlb_length) : n_regimes
    if get_setting(m, :gensys2)
        gensys_regimes = [1:first_gensys2_regime-1]
        if last_gensys2_regime != n_regimes
            append!(gensys_regimes, [last_gensys2_regime+1:n_regimes])
        end
    else
        gensys_regimes = [1:n_regimes]
    end
    gensys2_regimes = [first_gensys2_regime-1:last_gensys2_regime]

    # Solve model
    transitions = Vector{Vector{Transition{T}}}(undef, n_tvis)
    TTTs_vec    = Vector{Vector{Matrix{T}}}(undef, n_tvis)
    RRRs_vec    = Vector{Vector{Matrix{T}}}(undef, n_tvis)
    CCCs_vec    = Vector{Vector{Vector{T}}}(undef, n_tvis)
    for (i, replace_eqcond_func_dict) in enumerate(tvis_replace_eqcond_func_dict) # For each set of equilibrium conditions,
        if apply_altpolicy
            m <= Setting(:replace_eqcond_func_dict, replace_eqcond_func_dict)         # calculate the implied regime-switching system
        end

        TTTs_vec[i], RRRs_vec[i], CCCs_vec[i] = solve(m; apply_altpolicy = apply_altpolicy, regime_switching = regime_switching,
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
