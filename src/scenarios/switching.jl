"""
```
simulate_switching(m, scen::SwitchingScenario; verbose = :low)
```

Simulate switching in and out of a default scenario for the
SwitchingScenario scen. See `SwitchingScenario` for more
info. Returns a dictionary of results and the proportion
of times a switch actually occured.
"""
function simulate_switching(m::AbstractModel, scen::SwitchingScenario;
                            verbose::Symbol = :low)

    info_print(verbose, :low, "Simulating switching for " * string(scen.key) * "...")
    println(verbose, :low, "Start time: " * string(now()))
    println(verbose, :low, "Outputs will be saved in " * rawpath(m, "scenarios"))

    tic = time_ns()
    # Revert model alt policy to historical rule
    m <= Setting(:alternative_policy, AltPolicy(:historical, solve, eqcond), false, "apol",
                 "Alternative policy")

    # Get output file names
    original_output_files = get_scenario_output_files(m, scen.original,
                                                      [:forecastobs, :forecastpseudo])
    default_output_files = get_scenario_output_files(m, scen.default,
                                                     [:forecastobs, :forecastpseudo])

    results = Dict{Symbol, Array{Float64}}()
    switching_results = [0.0, 0.0]

    for (i, output_var) in enumerate([:forecastobs, :forecastpseudo])
        # Read in original and default draws
        original_draws = JLD2.jldopen(original_output_files[output_var], "r") do file
            read(file, "arr")
        end
        default_draws  = JLD2.jldopen(default_output_files[output_var], "r") do file
            read(file, "arr")
        end

        n_draws, n_vars, n_periods = size(original_draws)
        n_default_draws = size(default_draws, 1)

        @assert n_vars == size(default_draws, 2)
        @assert length(scen.probs_enter) == n_periods <= size(default_draws, 3)

        # Simulate switching
        results[output_var] = zeros(n_draws, n_vars, n_periods)
        proportion_switched = 0.0
        for j = 1:n_draws
            k = rand(1:n_default_draws)
            results[output_var][j, :, :], switched = switch(original_draws[j, :, :], default_draws[k, :, :],
                                                            scen.probs_enter, scen.probs_exit)
            proportion_switched += Float64(switched)
        end
        switching_results[i] = proportion_switched / n_draws
    end

    results[:proportion_switched] = [mean(switching_results)]

    # Write outputs
    output_files = get_scenario_output_files(m, scen, [:forecastobs, :forecastpseudo])
    write_scenario_forecasts(m, output_files, results, verbose = verbose)

    # Print
    switching_time = (time_ns() - tic)/1e9
    switching_time_min = switching_time/60
    println(verbose, :low, "\nTime elapsed: " * string(switching_time_min) * " minutes")
    println(verbose, :low, "Switching complete: " * string(now()))

    return results, switching_results
end

#for switching with multiple models
function simulate_switching(m_original::AbstractModel, m_default::AbstractModel, scen::SwitchingScenario;
                            verbose::Symbol = :low)

    info(verbose, :low, "Simulating switching for " * string(scen.key) * "...")
    println(verbose, :low, "Start time: " * string(now()))
    println(verbose, :low, "Outputs will be saved in " * rawpath(m_original, "scenarios"))
    tic()

    # Revert model alt policy to historical rule
    m_original <= Setting(:alternative_policy, AltPolicy(:historical, solve, eqcond), false, "apol",
                 "Alternative policy")
    m_default <= Setting(:alternative_policy, AltPolicy(:historical, solve, eqcond), false, "apol",
                 "Alternative policy")


    # Get output file names
    original_output_files = get_scenario_output_files(m_original, scen.original,
                                                      [:forecastobs, :forecastpseudo])
    default_output_files = get_scenario_output_files(m_default, scen.default,
                                                     [:forecastobs, :forecastpseudo])

    results = Dict{Symbol, Array{Float64}}()
    switching_results = [0.0, 0.0]

    for (i, output_var) in enumerate([:forecastobs, :forecastpseudo])
        if output_var==:forecastobs
            include_obs = intersect(keys(m_default.observables), keys(m_original.observables))
            inds_default = Vector(0)
            inds_original = Vector(0)
            for key in include_obs
                inds_default = vcat(inds_default, m_default.observables[key])
                inds_original = vcat(inds_original, m_original.observables[key])
            end
        elseif output_var==:forecastpseudo
            include_obs = intersect(keys(m_default.pseudo_observables), keys(m_original.pseudo_observables))
            inds_default = Vector(0)
            inds_original = Vector(0)
            for key in include_obs
                inds_default = vcat(inds_default, m_default.pseudo_observables[key])
                inds_original = vcat(inds_original, m_original.pseudo_observables[key])
            end
        end

        # Read in original and default draws
        original_draws = load(original_output_files[output_var], "arr")[:, inds_original, :]
        default_draws  = load(default_output_files[output_var], "arr")[:, inds_default, :]

        n_draws, n_vars, n_periods = size(original_draws)
        n_default_draws = size(default_draws, 1)

        @assert n_vars == size(default_draws, 2)
        @assert length(scen.probs_enter) == n_periods <= size(default_draws, 3)

        # Simulate switching
        results[output_var] = zeros(n_draws, n_vars, n_periods)
        proportion_switched = 0.0
        for j = 1:n_draws
            k = rand(1:n_default_draws)
            results[output_var][j, :, :], switched = switch(original_draws[j, :, :], default_draws[k, :, :],
                                                            scen.probs_enter, scen.probs_exit)
            proportion_switched += Float64(switched)
        end
        switching_results[i] = proportion_switched / n_draws
    end

    results[:proportion_switched] = [mean(switching_results)]

    # Write outputs
    output_files = get_scenario_output_files(m_original, scen, [:forecastobs, :forecastpseudo])
    write_scenario_forecasts(m_original, output_files, results, verbose = verbose)

    # Print
    switching_time = toq()
    switching_time_min = switching_time/60
    println(verbose, :low, "\nTime elapsed: " * string(switching_time_min) * " minutes")
    println(verbose, :low, "Switching complete: " * string(now()))

    return results, switching_results
end


"""
```
switch(original, default, probs_enter, probs_exit)
```

Simulate entry and exit from the original vector into the default one
according to the specified probabilities.
"""
function switch(original::Matrix{Float64}, default::Matrix{Float64},
                probs_enter::Vector{Float64}, probs_exit::Vector{Float64})

    n_vars, n_periods = size(original)
    output = zeros(n_vars, n_periods)

    # Find period of switch from default to original
    period_in = choose_switch_period(probs_enter, 1)

    # Record whether original is actually switched into
    switched = period_in < n_periods

    # Find period of switch from original to default
    # Start in period_in + 1 s.t. at least 1 period is spent in the original scenario
    period_out = choose_switch_period(probs_exit, period_in + 1)

    output[:, 1:(period_in-1)]          = default[:, 1:(period_in-1)]
    output[:, period_in:(period_out-1)] = original[:, 1:(period_out-period_in)]
    output[:, period_out:n_periods]     = default[:, period_out:n_periods]

    return output, switched
end

function choose_switch_period(probs::Vector{Float64}, first_period::Int)
    n_periods = length(probs)

    for t = first_period:n_periods
        if rand() < probs[t]
            return t
        end
    end

    return n_periods + 1
end
