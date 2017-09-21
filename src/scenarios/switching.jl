"""
```
simulate_switching(m, scen::SwitchingScenario; verbose = :low)
```

Simulate switching in and out of a default scenario for the
SwitchingScenario scen. See `SwitchingScenario` for more
info.
"""
function simulate_switching(m::AbstractModel, scen::SwitchingScenario;
                            verbose::Symbol = :low)


    if VERBOSITY[verbose] >= VERBOSITY[:low]
        info("Simulating switching for " * string(scen.key) * "...")
        println("Start time: " * string(now()))
        println("Outputs will be saved in " * rawpath(m, "scenarios"))
        tic()
    end

    original_output_files = get_scenario_output_files(m, scen.original,
                                                      [:forecastobs, :forecastpseudo])
    default_output_files = get_scenario_output_files(m, scen.default,
                                                     [:forecastobs, :forecastpseudo])

    results = Dict{Symbol, Array{Float64}}()
    for output_var in [:forecastobs, :forecastpseudo]
        # Read in original and default draws
        original_draws = load(original_output_files[output_var], "arr")
        default_draws  = load(default_output_files[output_var], "arr")

        n_draws, n_vars, n_periods = size(original_draws)
        n_default_draws = size(default_draws, 1)

        @assert n_vars == size(default_draws, 2)
        @assert length(scen.probs_enter) == n_periods <= size(default_draws, 3)

        # Simulate switching
        results[output_var] = zeros(n_draws, n_vars, n_periods)
        for i = 1:n_draws
            j = rand(1:n_default_draws)
            results[output_var][i, :, :] = switch(original_draws[i, :, :], default_draws[j, :, :],
                                                  scen.probs_enter, scen.probs_exit)
        end

    end

    # Write outputs
    output_files = get_scenario_output_files(m, scen, [:forecastobs, :forecastpseudo])
    write_scenario_forecasts(m, output_files, results, verbose = verbose)

    # Print
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        switching_time = toq()
        switching_time_min = switching_time/60
        println("\nTime elapsed: " * string(switching_time_min) * " minutes")
        println("Switching complete: " * string(now()))
    end

    return results
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

    # Find period of switch from original to default
    # Start in period_in + 1 s.t. at least 1 period is spent in the original scenario
    period_out = choose_switch_period(probs_exit, period_in + 1)

    output[:, 1:(period_in-1)]          = default[:, 1:(period_in-1)]
    output[:, period_in:(period_out-1)] = original[:, 1:(period_out-period_in)]
    output[:, period_out:n_periods]     = default[:, period_out:n_periods]

    return output
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