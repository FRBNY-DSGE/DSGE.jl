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


    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
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
        # these forecasts are n_draws x n_vars x n_periods
        original_draws = load(original_output_files[output_var])["arr"]
        default_draws  = load(default_output_files[output_var])["arr"]

        n_draws, n_vars, n_periods = size(original_draws)
        n_default_draws = size(default_draws, 1)

        @assert n_vars == size(default_draws, 2)
        @assert n_periods <= size(default_draws, 3)

        results[output_var] = zeros(n_draws, n_vars, n_periods)

        for i = 1:n_draws
            j = rand(1:n_default_draws)
            results[output_var][i, :, :] = switch(original_draws[i,:,:], default_draws[j,:,:],
                                                  scen.probs_enter, scen.probs_exit)
        end

    end

    output_files = get_scenario_output_files(m, scen, [:forecastobs, :forecastpseudo])
    write_scenario_forecasts(m, output_files, results, verbose = verbose)

    # Print
    if DSGE.VERBOSITY[verbose] >= DSGE.VERBOSITY[:low]
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

    # find period we switch into the original from the default
    period_in  = choose_last_period(probs_enter, 1, n_periods)
    # spend 1 period in the original scenario at least
    period_out = choose_last_period(probs_exit, period_in + 2, n_periods)

    output[:,1:period_in]            = default[:,1:period_in]
    output[:,period_in+1:period_out] = original[:,1:(period_out-period_in)]
    if period_out < n_periods
        output[:, period_out+1:n_periods] = default[:,period_out+1:n_periods]
    end

    return output
end

function choose_last_period(probs::Vector{Float64}, first_period::Int, n_periods::Int)
    last_period = n_periods

    for t = first_period:n_periods
        if rand() < probs[t]
            last_period = t-1
            return last_period
        end
    end

    return last_period
end