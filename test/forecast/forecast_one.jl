using DSGE, Base.Test, JLD, DistributedArrays
include("../util.jl")

path = dirname(@__FILE__)

# Initialize model object
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks    => Setting(:n_anticipated_shocks, 6),
    :use_population_forecast => Setting(:use_population_forecast, true),
    :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :date_conditional_end    => Setting(:date_conditional_end, quartertodate("2015-Q4")),
    :date_forecast_end       => Setting(:date_forecast_end, quartertodate("2016-Q1")),
    :forecast_kill_shocks    => Setting(:forecast_kill_shocks, true),
    :saveroot                => Setting(:saveroot, normpath(joinpath(dirname(@__FILE__), "..", "reference"))),
    :use_parallel_workers    => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

# Add parallel workers
my_procs = addprocs(5)
@everywhere using DSGE

# Run forecasts
forecast_outputs = Dict{Tuple{Symbol, Symbol}, Dict{Symbol, Any}}()
output_vars = [:histstates, :histpseudo, :histshocks,
               :forecaststates, :forecastpseudo, :forecastobs, :forecastshocks,
               :shockdecstates, :shockdecpseudo, :shockdecobs]
output_files = []

# Call forecast_one once without timing
df = load_data(m; verbose = :none)
DSGE.compile_forecast_one(m, df; cond_type = :none, output_vars = output_vars,
             verbose = :none, procs = my_procs)

# Check error handling for input_type = :subset
@test_throws ErrorException forecast_one(m, df; input_type = :subset, cond_type = :none,
                                output_vars = output_vars, subset_inds = collect(1:10),
                                subset_string = "", verbose = :none)
@test_throws ErrorException forecast_one(m, df; input_type = :subset, cond_type = :none,
                                output_vars = output_vars, subset_inds = Vector{Int}(),
                                subset_string = "test", verbose = :none)

# Run all forecast combinations
for input_type in [:init, :mode, :full]

    df = load_data(m; verbose = :none)

    for cond_type in [:none, :semi, :full]

        println("input_type = $(input_type), cond_type = $(cond_type)")

        forecast_output = Dict{Symbol, Any}()
        forecast_output[:df] = load_data(m; cond_type=cond_type, try_disk=true, verbose=:none)

        procs = input_type == :full ? my_procs : [myid()]
        @time new_forecast = forecast_one(m, forecast_output[:df];
            input_type = input_type, cond_type = cond_type, output_vars = output_vars,
            verbose = :none, procs = procs)
        merge!(forecast_output, new_forecast)

        forecast_outputs[(cond_type, input_type)] = forecast_output

        new_output_files = collect(values(DSGE.get_output_files(m, input_type, output_vars, cond_type)))
        append!(output_files, new_output_files)
    end
end

# Read expected output
exp_forecast_outputs = jldopen("$path/../reference/forecast_one_out.jld", "r") do file
    read(file, "exp_forecast_outputs")
end

# Test forecast outputs
for input_type in [:init, :mode]

    if input_type == :init
        DSGE.init_parameters!(m)
        DSGE.steadystate!(m)
    elseif input_type == :mode
        specify_mode!(m, DSGE.get_input_file(m, :mode))
    end

    for cond_type in [:none, :semi, :full]

        ## Historical states and shocks
        histstates = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:histstates], 1, :, :))
        histshocks = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:histshocks], 1, :, :))
        histpseudo = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:histpseudo], 1, :, :))

        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:histstates] histstates
        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:histshocks] histshocks
        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:histpseudo] histpseudo

        ## Forecasted states, observables, pseudos, and shocks
        forecaststates = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:forecaststates], 1, :, :))
        forecastobs    = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:forecastobs], 1, :, :))
        forecastpseudo = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:forecastpseudo], 1, :, :))
        forecastshocks = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:forecastshocks], 1, :, :))

        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:forecaststates] forecaststates
        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:forecastobs]    forecastobs
        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:forecastpseudo] forecastpseudo
        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:forecastshocks] forecastshocks

        ## Shock decompositions of states, pseudo-observables, and observables
        shockdecstates = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:shockdecstates], 1, :, :, :))
        shockdecobs    = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:shockdecobs], 1, :, :, :))
        shockdecpseudo = convert(Array, slice(forecast_outputs[(cond_type, input_type)][:shockdecpseudo], 1, :, :, :))

        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:shockdecstates] shockdecstates
        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:shockdecobs]    shockdecobs
        @test_matrix_approx_eq exp_forecast_outputs[(cond_type, input_type)][:shockdecpseudo] shockdecpseudo

    end # cond_type

end # input_type


# Delete all files written by forecast_one
map(rm, unique(output_files))

# Remove parallel workers
rmprocs(my_procs)

nothing
