include("drivers.jl")

"""
```
forecast_one_new(m::AbstractModel, df::DataFrame; input_type::Symbol = :mode,
    output_type::Symbol = :simple, cond_type::Symbol = :none)
```

New version of `forecast_one`, created for the purposes of benchmarking
performance against the existing `forecast_one`. Eventually these changes will
be merged back into `forecast_one`.
"""
function forecast_one_new(m::AbstractModel, df::DataFrame;
                          input_type::Symbol = :mode,
                          output_vars::Vector{Symbol} = [],
                          cond_type::Symbol = :none,
                          verbose::Symbol = :low)

    ### 1. Setup

    # Prepare forecast inputs
    systems, states = prepare_forecast_inputs(m, df; input_type = input_type,
        cond_type = cond_type)
    ndraws = length(systems)

    # Prepare forecast outputs
    forecast_output = Dict{Symbol, Array{Float64}}()
    forecast_output_files = get_output_files(m, input_type, output_vars, cond_type)
    output_dir = rawpath(m, "forecast")
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nForecasting input_type = $input_type, cond_type = $cond_type...")
        println("Forecast outputs will be saved in $output_dir")
    end

    # Inline definition s.t. the dicts forecast_output and forecast_output_files are accessible
    function write_forecast_outputs(vars::Vector{Symbol})
        for var in vars
            file = forecast_output_files[var]
            jldopen(file, "w") do f
                write(f, string(var), forecast_output[var])
                if VERBOSITY[verbose] >= VERBOSITY[:high]
                    println(" * Wrote $(basename(file))")
                end
            end
        end
    end


    ### 2. Smoothed Histories

    # must re-run filter/smoother for conditional data in addition to explicit cases
    hist_vars = [:histstates, :histpseudo, :histshocks]
    shockdec_vars = [:shockdecstates, :shockdecpseudo, :shockdecobs]
    filterandsmooth_vars = vcat(hist_vars, shockdec_vars)

    if !isempty(intersect(output_vars, filterandsmooth_vars)) || cond_type in [:semi, :full]

        histstates, histshocks, histpseudo, kals = filterandsmooth(m, df, systems; cond_type = cond_type)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))

            forecast_output[:histstates] = histstates[:, 1:T, :]
            forecast_output[:histshocks] = histshocks[:, 1:T, :]
            forecast_output[:histpseudo] = histpseudo[:, 1:T, :]
        else
            forecast_output[:histstates] = histstates
            forecast_output[:histshocks] = histshocks
            forecast_output[:histpseudo] = histpseudo
        end

        write_forecast_outputs(hist_vars)
    end


    ### 3. Forecasts

    # For conditional data, use the end of the hist states as the initial state
    # vector for the forecast
    if cond_type in [:semi, :full]
        states = [kal[:zend]::Vector{Float64} for kal in kals]
    end

    forecast_vars = [:forecaststates, :forecastobs, :forecastpseudo, :forecastshocks]

    if !isempty(intersect(output_vars, forecast_vars))
        forecaststates, forecastobs, forecastpseudo, forecastshocks =
            forecast(m, systems, states)

        # For conditional data, transplant the obs/state/pseudo vectors from hist to forecast
        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))
            histobs = df_to_matrix(m, df; cond_type = cond_type)[:, index_prezlb_start(m):end]
            histobs = repeat(histobs, outer = [1, 1, ndraws])

            forecast_output[:forecaststates] = cat(2, histstates[:, T+1:end, :], forecaststates)
            forecast_output[:forecastshocks] = cat(2, histshocks[:, T+1:end, :], forecastshocks)
            forecast_output[:forecastpseudo] = cat(2, histpseudo[:, T+1:end, :], forecastpseudo)
            forecast_output[:forecastobs]    = cat(2, histobs[:,    T+1:end, :], forecastobs)
        else
            forecast_output[:forecaststates] = forecaststates
            forecast_output[:forecastshocks] = forecastshocks
            forecast_output[:forecastpseudo] = forecastpseudo
            forecast_output[:forecastobs]    = forecastobs
        end

        write_forecast_outputs(forecast_vars)
    end


    ### 4. Shock Decompositions

    if !isempty(intersect(output_vars, [:shockdecstates, :shockdecobs, :shockdecpseudo]))
        histshocks = [histshocks[:, :, i]::Matrix{Float64} for i = 1:ndraws]
        shockdecstates, shockdecobs, shockdecpseudo = shock_decompositions(m, systems, histshocks)

        forecast_output[:shockdecstates] = shockdecstates
        forecast_output[:shockdecpseudo] = shockdecpseudo
        forecast_output[:shockdecobs]    = shockdecobs

        write_forecast_outputs(shockdec_vars)
    end


    # Return only saved elements of dict
    filter!((k, v) -> k ∈ output_vars, forecast_output)
    return forecast_output
end
