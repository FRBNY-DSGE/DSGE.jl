"""
```
compute_meansbands(m, input_type, cond_type, output_vars; forecast_string = "",
    verbose = :low, kwargs...)

compute_meansbands(m, input_type, cond_type, output_var, df; forecast_string = "",
    population_data = DataFrame(), population_forecast = DataFrame(),
    verbose = :none, kwargs...)

compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
    forecast_string = "", population_data = DataFrame(),
    population_forecast = DataFrame(), verbose = :low,
    kwargs...)
```

Compute means and bands for pseudo-observables, observables, and shocks, and write
the results to a file. Other methods are for one `output_var` and one `var_name` respectively.

### Keyword Arguments

- `forecast_string::String`: forecast identifier string (the value
  \"fcid=value\" in the forecast output filename). Required when
  `input_type == :subset`

- `density_bands::Vector{Float64}`: a vector of percent values (between 0 and 1) for
  which to compute density bands

- `minimize::Bool`: if `true`, choose shortest interval, otherwise just chop off
  lowest and highest (percent/2)

- `verbose`: level of error messages to be printed to screen. One of `:none`,
  `:low`, `:high`

- `check_empty_columns::Bool`: if true, throw an error if
    calling load_data yields an empty column.
"""
function compute_meansbands(m::AbstractDSGEModel, input_type::Symbol,
                            cond_type::Symbol, output_vars::Vector{Symbol};
                            forecast_string::String = "",
                            verbose::Symbol = :low, df::DataFrame = DataFrame(),
                            check_empty_columns::Bool = true,
                            bdd_fcast::Bool = true,
                            kwargs...)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        output_dir = workpath(m, "forecast")
        println()
        @Base.info "Computing means and bands for input_type = $input_type, cond_type = $cond_type..."
        println("Start time: $(now())")
        println("Means and bands will be saved in $output_dir")
    end
    elapsed_time = @elapsed let
        # Determine full set of output_vars necessary for plotting desired result
        output_vars = add_requisite_output_vars(output_vars; bdd_fcast = bdd_fcast)
        if input_type == :prior
            output_vars = setdiff(output_vars, [:bddforecastobs])
        end
        # Load population data and main dataset (required for some transformations)
        if all(var -> get_product(var) == :irf, output_vars)
            population_data, population_forecast = DataFrame(), DataFrame()
        else
            population_data, population_forecast = load_population_growth(m, verbose = verbose)
            isempty(df) && (df = load_data(m, cond_type = cond_type,
                                           check_empty_columns = check_empty_columns, verbose = :none))
        end
        for output_var in output_vars
            prod = get_product(output_var)
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                if prod in [:shockdec, :irf]
                    println("Computing " * string(output_var) * " for shocks:")
                else
                    print("Computing " * string(output_var) * "... ")
                end
            end

            # Compute means and bands
            mb = compute_meansbands(m, input_type, cond_type, output_var, df;
                                    forecast_string = forecast_string,
                                    population_data = population_data,
                                    population_forecast = population_forecast,
                                    verbose = verbose,
                                    kwargs...)
            GC.gc()
        end
    end
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        total_mb_time     = elapsed_time
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function compute_meansbands(m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                            output_var::Symbol, df::DataFrame;
                            forecast_string::String = "",
                            population_data::DataFrame = DataFrame(),
                            population_forecast::DataFrame = DataFrame(),
                            verbose::Symbol = :none,
                            kwargs...)

    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast metadata
    metadata = get_mb_metadata(m, input_type, cond_type, output_var; forecast_string = forecast_string)

    date_list      = product == :irf ? Date[] : collect(keys(metadata[:date_inds]))
    variable_names = collect(keys(metadata[:indices]))
    pop_growth     = get_mb_population_series(product, population_data, population_forecast, date_list)

    # Compute means and bands
    if product in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                   :bddforecast, :bddforecastut, :bddforecast4q, :dettrend, :trend]
                   # :forecastlvl, :histlvl, :bddforecastlvl]
        # Get to work!
        # pmap produces an error for trendobs sometimes, so just doing this iteratively
        if product == :trend
            mb_vec = Vector{Any}(undef,length(variable_names))
            for i in 1:length(mb_vec)
                mb_vec[i] = compute_meansbands(m, input_type, cond_type, output_var,
                                               variable_names[i], df; pop_growth = pop_growth,
                                               forecast_string = forecast_string,
                                               kwargs...)
            end
        else
            mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
                                                         pop_growth = pop_growth, forecast_string = forecast_string, kwargs...),
                          variable_names)
        end

        # Re-assemble pmap outputs
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
            means[!,var_name] = var_means
            bands[var_name] = var_bands
            bands[var_name][!,:date] = date_list
        end

    elseif product in [:shockdec, :irf]
        means = product == :irf ? DataFrame() : DataFrame(date = date_list)
        bands = Dict{Symbol, DataFrame}()

        # Get to work!
        for shock_name in keys(metadata[:shock_indices])
            println(verbose, :high, "  * " * string(shock_name))

            mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
                                          pop_growth = pop_growth, shock_name = Nullables.Nullable(shock_name),
                                                         forecast_string = forecast_string,
                                                         kwargs...),
                          variable_names)

            # Re-assemble pmap outputs
            for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
                means[!, Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_means
                bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_bands
                if product != :irf
                    bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)][!, :date] = date_list
                end
            end
        end

    else
        error("Invalid product: $product")
    end

    mb = MeansBands(metadata, means, bands)

    # Write to file
    filepath = get_meansbands_output_file(m, input_type, cond_type, output_var,
                                          forecast_string = forecast_string)
    dirpath = dirname(filepath)
    isdir(dirpath) || mkpath(dirpath)
    JLD2.jldopen(filepath, true, true, true, IOStream) do file
        write(file, "mb", mb)
    end

    sep = prod in [:shockdec, :irf] ? "  " : ""
    println(verbose, :high, sep * "wrote " * basename(filepath))

    return mb
end

function compute_meansbands(m::AbstractDSGEModel, input_type::Symbol, cond_type::Symbol,
                            output_var::Symbol, var_name::Symbol, df::DataFrame;
                            forecast_string::String = "",
                            pop_growth::AbstractVector{Float64} = Float64[],
                            shock_name::Nullable{Symbol} = Nullables.Nullable{Symbol}(),
                            density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                            minimize::Bool = false,
                            compute_shockdec_bands::Bool = false)

    # Return only one set of bands if we read in only one draw
    if input_type in [:init, :mode, :mean]
        density_bands = [.5]
    end

    # Determine class and product
    class = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast draws
    fcast_series, transform = read_forecast_output(m, input_type, cond_type,
                                                   output_var, var_name, shock_name,
                                                   forecast_string = forecast_string)

    # Reverse transform
    y0_index = get_y0_index(m, product) # this should be index_forecast_start(m) - 4, so 4 quarters before forecast
    data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,var_name]), NaN))) : fill(NaN, size(df, 1))
    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))
    transformed_series = mb_reverse_transform(fcast_series, transform, product, class,
                                              y0_index = y0_index, data = data,
                                              pop_growth = pop_growth)

    # Handle NaNs
    if output_var != :histobs && any(isnan.(transformed_series))
        # Remove rows with NaNs
        nanrows = vec(mapslices(x -> any(isnan.(x)), transformed_series, dims = Int[2]))
        transformed_series = transformed_series[.!nanrows, :]
    end

    # Compute means and bands
    means = vec(mean(transformed_series, dims = 1))
    bands = if product in [:shockdec, :dettrend, :trend] && !compute_shockdec_bands
        Dict{Symbol,DataFrame}()
    else
        find_density_bands(transformed_series, density_bands, minimize = minimize)
    end
    return means, bands
end

function compute_meansbands(models::Vector,
                            input_types::Vector{Symbol},
                            cond_types::Vector{Symbol},
                            output_vars::Vector{Symbol};
                            forecast_strings::Vector{String} = Vector{String}("", 0),
                            weights::Vector{Float64} = Vector{Float64}(undef, 0),
                            combo_forecast_string::String = join(forecast_strings) * join(map(x->string(x), weights), "_"),
                            df::DataFrame = DataFrame(),
                            check_empty_columns::Bool = true,
                            variable_names::Vector{Symbol} = Vector{Symbol}(undef, 0),
                            shock_names::Vector{Symbol} = Vector{Symbol}(undef, 0),
                            verbose::Symbol = :low,
                            kwargs...)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        output_dir = workpath(models[1], "forecast")
        println()
        @Base.info "Computing means and bands for input_types = $input_types, cond_types = $cond_types..."
        println("Start time: $(now())")
        println("Means and bands will be saved in $output_dir")
    end
    elapsed_time = @elapsed let
        # Determine full set of output_vars necessary for plotting desired result
        output_vars = add_requisite_output_vars(output_vars)
        if input_types[1] == :prior
            output_vars = setdiff(output_vars, [:bddforecastobs])
        end
        # Load population data and main dataset (required for some transformations)
        if all(var -> get_product(var) == :irf, output_vars)
            population_data, population_forecast = DataFrame(), DataFrame()
        else
            population_data, population_forecast = load_population_growth(models[1], verbose = verbose)
            isempty(df) && (df = load_data(models[1], check_empty_columns = check_empty_columns, verbose = :none))
        end
        for output_var in output_vars
            prod = get_product(output_var)
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                if prod in [:shockdec, :irf]
                    println("Computing " * string(output_var) * " for shocks:")
                else
                    print("Computing " * string(output_var) * "... ")
                end
            end

            # Compute means and bands
            mb = compute_meansbands(models, input_types,
                                    cond_types, output_var, df;
                                    forecast_strings = forecast_strings,
                                    combo_forecast_string = combo_forecast_string,
                                    weights = weights,
                                    population_data = population_data,
                                    population_forecast = population_forecast,
                                    variable_names = variable_names,
                                    shock_names = shock_names,
                                    verbose = verbose,
                                    kwargs...)
            GC.gc()
        end
    end
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        total_mb_time     = elapsed_time
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function compute_meansbands(models::Vector,
                            input_types::Vector{Symbol},
                            cond_types::Vector{Symbol},
                            output_var::Symbol, df::DataFrame;
                            forecast_strings::Vector{String} = Vector{String}("", 0),
                            combo_forecast_string::String = "",
                            weights::Vector{Float64} = Vector{Float64}(undef, 0),
                            population_data::DataFrame = DataFrame(),
                            population_forecast::DataFrame = DataFrame(),
                            variable_names::Vector{Symbol} = Vector{Symbol}(undef, 0),
                            shock_names::Vector{Symbol} = Vector{Symbol}(undef, 0),
                            verbose::Symbol = :none,
                            kwargs...)

    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast metadata
    metadata = get_mb_metadata(models[1], input_types[1], cond_types[1], output_var; forecast_string = forecast_strings[1])
    metadata[:forecast_string] = combo_forecast_string

    date_list      = product == :irf ? Date[] : collect(keys(metadata[:date_inds]))
    if isempty(variable_names)
        variable_names = collect(keys(metadata[:indices]))
    end
    pop_growth     = get_mb_population_series(product, population_data, population_forecast, date_list)

    # Compute means and bands
    if product in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                   :bddforecast, :bddforecastut, :bddforecast4q, :dettrend, :trend]
        # Get to work!
        # pmap produces an error for trendobs sometimes, so just doing this iteratively
        if product == :trend
            mb_vec = Vector{Any}(undef,length(variable_names))
            for i in 1:length(mb_vec)
                mb_vec[i] = compute_meansbands(models, input_types,
                                               cond_types, output_var,
                                               variable_names[i], df; pop_growth = pop_growth,
                                               forecast_strings = forecast_strings,
                                               weights = weights,
                                               kwargs...)
            end
        else
            mb_vec = pmap(var_name -> compute_meansbands(models, input_types, cond_types, output_var, var_name, df;
                                                         pop_growth = pop_growth,
                                                         forecast_strings = forecast_strings, weights = weights, kwargs...),
                          variable_names)
        end

        # Re-assemble pmap outputs
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
            means[!,var_name] = vec(var_means)
            bands[var_name] = var_bands
            bands[var_name][!,:date] = date_list
        end

    elseif product in [:shockdec, :irf]
        means = product == :irf ? DataFrame() : DataFrame(date = date_list)
        bands = Dict{Symbol, DataFrame}()

        # Get to work!
        if isempty(shock_names)
            shock_names = keys(metadata[:shock_indices])
        end
        for shock_name in shock_names
            println(verbose, :high, "  * " * string(shock_name))

            mb_vec = pmap(var_name -> compute_meansbands(models, input_types,
                                           cond_types, output_var,
                                           var_name, df; pop_growth = pop_growth,
                                           forecast_strings = forecast_strings,
                                           weights = weights, shock_name = Nullables.Nullable(shock_name),
                                           kwargs...),
                             variable_names)

            # Re-assemble pmap outputs
            for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
                means[!, Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_means
                bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_bands
                if product != :irf
                    bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)][!, :date] = date_list
                end
            end
        end

    else
        error("Invalid product: $product")
    end # of if product

    mb = MeansBands(metadata, means, bands)

    # Write to file
    filepath = get_meansbands_output_file(models[1], input_types[1], cond_types[1], output_var,
                                          forecast_string = combo_forecast_string)

    dirpath = dirname(filepath)
    isdir(dirpath) || mkpath(dirpath)
    JLD2.jldopen(filepath, true, true, true, IOStream) do file
        write(file, "mb", mb)
    end

    sep = prod in [:shockdec, :irf] ? "  " : ""
    println(verbose, :high, sep * "wrote " * basename(filepath))

    return mb
end

function compute_meansbands(models::Vector,
                            input_types::Vector{Symbol},
                            cond_types::Vector{Symbol},
                            output_var::Symbol, var_name::Symbol, df::DataFrame;
                            forecast_strings::Vector{String} = Vector{String}("", 0),
                            weights::Vector{Float64} = Vector{Float64}(undef, 0),
                            pop_growth::AbstractVector{Float64} = Float64[],
                            shock_name::Nullable{Symbol} = Nullables.Nullable{Symbol}(),
                            density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                            minimize::Bool = false,
                            compute_shockdec_bands::Bool = false)

    # Return only one set of bands if we read in only one draw
    if input_types[1] in [:init, :mode, :mean]
        density_bands = [.5]
    end

    # Determine class and product
    class = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast draws
    fcast_seriess = Vector{Matrix{Float64}}(undef, length(models))
    transforms = Vector(undef, length(models))
    for i in 1:length(models)
        fcast_seriess[i], transforms[i] = read_forecast_output(models[i], input_types[i], cond_types[i],
                                                         output_var, var_name, shock_name,
                                                         forecast_string = forecast_strings[i])
    end

    # Reverse transform
    y0_index = get_y0_index(models[1], product)
    data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,var_name]), NaN))) : fill(NaN, size(df, 1))
    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))

    fcast_series = Matrix{Float64}(undef, 0, size(fcast_seriess[1], 2))
    if all(map(x->size(x, 1) > 1, fcast_seriess))
        for i in 1:length(models)
            sampled = sample(1:size(fcast_seriess[i], 1), Int(20000*weights[i]))
            fcast_series = vcat(fcast_series, fcast_seriess[i][sampled, :])
        end
    else
        fcast_series = zeros(size(fcast_seriess[1]))
        for i in 1:length(models)
            fcast_series = fcast_series + weights[i]*fcast_seriess[i]
        end
    end

    transformed_series = mb_reverse_transform(fcast_series, transforms[1], product, class,
                                              y0_index = y0_index, data = data,
                                              pop_growth = pop_growth)

    # Compute means and bands
    means = vec(mean(transformed_series, dims= 1))
    bands = if product in [:shockdec, :dettrend, :trend] && !compute_shockdec_bands
        Dict{Symbol,DataFrame}()
    else
        find_density_bands(transformed_series, density_bands, minimize = minimize)
    end
    return means, bands
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type1::Symbol, input_type2::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_vars::Vector{Symbol};
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            verbose::Symbol = :low, df::DataFrame = DataFrame(),
                            check_empty_columns::Bool = true,
                            kwargs...)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        output_dir = workpath(m1, "forecast")
        println()
        @Base.info "Computing means and bands for input_type = $input_type, cond_type = $cond_type..."
        println("Start time: $(now())")
        println("Means and bands will be saved in $output_dir")
    end
    elapsed_time = @elapsed let
        # Determine full set of output_vars necessary for plotting desired result
        output_vars = add_requisite_output_vars(output_vars)
        if input_type1 == :prior
            output_vars = setdiff(output_vars, [:bddforecastobs])
        end
        # Load population data and main dataset (required for some transformations)
        if all(var -> get_product(var) == :irf, output_vars)
            population_data, population_forecast = DataFrame(), DataFrame()
        else
            population_data, population_forecast = load_population_growth(m1, verbose = verbose)
            isempty(df) && (df = load_data(m1, check_empty_columns = check_empty_columns, verbose = :none))
        end
        for output_var in output_vars
            prod = get_product(output_var)
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                if prod in [:shockdec, :irf]
                    println("Computing " * string(output_var) * " for shocks:")
                else
                    print("Computing " * string(output_var) * "... ")
                end
            end

            # Compute means and bands
            mb = compute_meansbands(m1, m2, input_type1, input_type2,
                                    cond_type1, cond_type2, output_var, df;
                                    forecast_string1 = forecast_string1,
                                    forecast_string2 = forecast_string2,
                                    population_data = population_data,
                                    population_forecast = population_forecast,
                                    verbose = verbose,
                                    kwargs...)
            GC.gc()
        end
    end
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        total_mb_time     = elapsed_time
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type1::Symbol, input_type2::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_var::Symbol, df::DataFrame;
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            population_data::DataFrame = DataFrame(),
                            population_forecast::DataFrame = DataFrame(),
                            verbose::Symbol = :none,
                            kwargs...)

    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast metadata
    metadata = get_mb_metadata(m1, input_type1, cond_type1, output_var; forecast_string = forecast_string1)

    date_list      = product == :irf ? Date[] : collect(keys(metadata[:date_inds]))
    variable_names = collect(keys(metadata[:indices]))
    pop_growth     = get_mb_population_series(product, population_data, population_forecast, date_list)

    # Compute means and bands
    if product in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                   :bddforecast, :bddforecastut, :bddforecast4q, :dettrend, :trend]
        # Get to work!
        # pmap produces an error for trendobs sometimes, so just doing this iteratively
        mb_vec = Vector{Any}(undef,length(variable_names))
        for i in 1:length(mb_vec)
            mb_vec[i] = compute_meansbands(m1, m2, input_type1, input_type2,
                                           cond_type1, cond_type2, output_var,
                                           variable_names[i], df; pop_growth = pop_growth,
                                           forecast_string1 = forecast_string1,
                                           forecast_string2 = forecast_string2,
                                           kwargs...)
        end
        # mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
        #                               pop_growth = pop_growth, forecast_string = forecast_string, kwargs...),
        #               variable_names)

        # Re-assemble pmap outputs
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
            means[!,var_name] = var_means
            bands[var_name] = var_bands
            bands[var_name][!,:date] = date_list
        end

    elseif product in [:shockdec, :irf]
        means = product == :irf ? DataFrame() : DataFrame(date = date_list)
        bands = Dict{Symbol, DataFrame}()

        # Get to work!
        for shock_name in keys(metadata[:shock_indices])
            println(verbose, :high, "  * " * string(shock_name))

            mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
                                          pop_growth = pop_growth, shock_name = Nullables.Nullable(shock_name),
                                          forecast_string1 = forecast_string1,
                                          forecast_string2 = forecast_string2,
                                                         kwargs...),
                          variable_names)

            # Re-assemble pmap outputs
            for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
                means[!, Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_means
                bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_bands
                if product != :irf
                    bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)][!, :date] = date_list
                end
            end
        end

    else
        error("Invalid product: $product")
    end # of if product

    mb = MeansBands(metadata, means, bands)

    # Write to file
    filepath = get_meansbands_output_file(m1, input_type1, cond_type1, output_var,
                                          forecast_string = forecast_string1*forecast_string2)
    dirpath = dirname(filepath)
    isdir(dirpath) || mkpath(dirpath)
    JLD2.jldopen(filepath, true, true, true, IOStream) do file
        write(file, "mb", mb)
    end

    sep = prod in [:shockdec, :irf] ? "  " : ""
    println(verbose, :high, sep * "wrote " * basename(filepath))

    return mb
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type1::Symbol, input_type2::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_var::Symbol, var_name::Symbol, df::DataFrame;
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            pop_growth::AbstractVector{Float64} = Float64[],
                            shock_name::Nullable{Symbol} = Nullables.Nullable{Symbol}(),
                            density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                            minimize::Bool = false,
                            compute_shockdec_bands::Bool = false)

    # Return only one set of bands if we read in only one draw
    if input_type1 in [:init, :mode, :mean]
        density_bands = [.5]
    end

    # Determine class and product
    class = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast draws
    fcast_series1, transform1 = read_forecast_output(m1, input_type1, cond_type1,
                                                   output_var, var_name, shock_name,
                                                   forecast_string = forecast_string1)
    fcast_series2, transform2 = read_forecast_output(m2, input_type2, cond_type2,
                                                   output_var, var_name, shock_name,
                                                   forecast_string = forecast_string2)

    # Reverse transform
    y0_index = get_y0_index(m1, product)
    yt_index = get_yt_index(m1, product)
    data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,var_name]), NaN))) : fill(NaN, size(df, 1))
    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))

   #= @show size(fcast_series1)
    @show size(fcast_series2)
    if size(fcast_series1, 1) < size(fcast_series2, 1)
        mult = Int(size(fcast_series2, 1) / size(fcast_series1, 1))
        fcast_series1 = repeat(fcast_series1, mult, 1)
    else
        mult = Int(size(fcast_series1, 1) / size(fcast_series2, 1))
        fcast_series2 = repeat(fcast_series2, mult, 1)
    end
    @show size(fcast_series1)
    @show size(fcast_series2)

    fcast_series = 0.25*fcast_series1 + 0.75*fcast_series2=#

    if size(fcast_series1, 1) > 1 && size(fcast_series2, 1) > 1
        sampled1 = sample(1:size(fcast_series1, 1), Int(20000*.1))
        sampled2 = sample(1:size(fcast_series2, 1), Int(20000*.9))
        fcast_series = vcat(fcast_series1[sampled1, :], fcast_series2[sampled2, :])
    else
        fcast_series = .1*fcast_series1 + .9*fcast_series2
    end


    transformed_series = mb_reverse_transform(fcast_series, transform1, product, class,
                                              y0_index = y0_index, yt_index = yt_index,
                                              data = data,
                                              pop_growth = pop_growth)

    # Compute means and bands
    means = vec(mean(transformed_series, dims= 1))
    bands = if product in [:shockdec, :dettrend, :trend] && !compute_shockdec_bands
        Dict{Symbol,DataFrame}()
    else
        find_density_bands(transformed_series, density_bands, minimize = minimize)
    end
    return means, bands
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type1::Symbol, input_type2::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_vars::Vector{Symbol};
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            verbose::Symbol = :low, df::DataFrame = DataFrame(),
                            check_empty_columns::Bool = true,
                            kwargs...)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        output_dir = workpath(m1, "forecast")
        println()
        @Base.info "Computing means and bands for input_type = $input_type, cond_type = $cond_type..."
        println("Start time: $(now())")
        println("Means and bands will be saved in $output_dir")
    end
    elapsed_time = @elapsed let
        # Determine full set of output_vars necessary for plotting desired result
        output_vars = add_requisite_output_vars(output_vars)
        if input_type1 == :prior
            output_vars = setdiff(output_vars, [:bddforecastobs])
        end
        # Load population data and main dataset (required for some transformations)
        if all(var -> get_product(var) == :irf, output_vars)
            population_data, population_forecast = DataFrame(), DataFrame()
        else
            population_data, population_forecast = load_population_growth(m1, verbose = verbose)
            isempty(df) && (df = load_data(m1, check_empty_columns = check_empty_columns, verbose = :none))
        end
        for output_var in output_vars
            prod = get_product(output_var)
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                if prod in [:shockdec, :irf]
                    println("Computing " * string(output_var) * " for shocks:")
                else
                    print("Computing " * string(output_var) * "... ")
                end
            end

            # Compute means and bands
            mb = compute_meansbands(m1, m2, input_type1, input_type2,
                                    cond_type1, cond_type2, output_var, df;
                                    forecast_string1 = forecast_string1,
                                    forecast_string2 = forecast_string2,
                                    population_data = population_data,
                                    population_forecast = population_forecast,
                                    verbose = verbose,
                                    kwargs...)
            GC.gc()
        end
    end
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        total_mb_time     = elapsed_time
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type1::Symbol, input_type2::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_var::Symbol, df::DataFrame;
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            population_data::DataFrame = DataFrame(),
                            population_forecast::DataFrame = DataFrame(),
                            verbose::Symbol = :none,
                            kwargs...)

    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast metadata
    metadata = get_mb_metadata(m1, input_type1, cond_type1, output_var; forecast_string = forecast_string1)

    date_list      = product == :irf ? Date[] : collect(keys(metadata[:date_inds]))
    variable_names = collect(keys(metadata[:indices]))
    pop_growth     = get_mb_population_series(product, population_data, population_forecast, date_list)

    # Compute means and bands
    if product in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                   :bddforecast, :bddforecastut, :bddforecast4q, :dettrend, :trend]
        # Get to work!
        # pmap produces an error for trendobs sometimes, so just doing this iteratively
        mb_vec = Vector{Any}(undef,length(variable_names))
        for i in 1:length(mb_vec)
            mb_vec[i] = compute_meansbands(m1, m2, input_type1, input_type2,
                                           cond_type1, cond_type2, output_var,
                                           variable_names[i], df; pop_growth = pop_growth,
                                           forecast_string1 = forecast_string1,
                                           forecast_string2 = forecast_string2,
                                           kwargs...)
        end
        # mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
        #                               pop_growth = pop_growth, forecast_string = forecast_string, kwargs...),
        #               variable_names)

        # Re-assemble pmap outputs
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
            means[!,var_name] = var_means
            bands[var_name] = var_bands
            bands[var_name][!,:date] = date_list
        end

    elseif product in [:shockdec, :irf]
        means = product == :irf ? DataFrame() : DataFrame(date = date_list)
        bands = Dict{Symbol, DataFrame}()

        # Get to work!
        for shock_name in keys(metadata[:shock_indices])
            println(verbose, :high, "  * " * string(shock_name))

            mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
                                          pop_growth = pop_growth, shock_name = Nullables.Nullable(shock_name),
                                          forecast_string1 = forecast_string1,
                                          forecast_string2 = forecast_string2,
                                                         kwargs...),
                          variable_names)

            # Re-assemble pmap outputs
            for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
                means[!, Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_means
                bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_bands
                if product != :irf
                    bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)][!, :date] = date_list
                end
            end
        end

    else
        error("Invalid product: $product")
    end # of if product

    mb = MeansBands(metadata, means, bands)

    # Write to file
    filepath = get_meansbands_output_file(m1, input_type1, cond_type1, output_var,
                                          forecast_string = forecast_string1*forecast_string2)
    dirpath = dirname(filepath)
    isdir(dirpath) || mkpath(dirpath)
    JLD2.jldopen(filepath, true, true, true, IOStream) do file
        write(file, "mb", mb)
    end

    sep = prod in [:shockdec, :irf] ? "  " : ""
    println(verbose, :high, sep * "wrote " * basename(filepath))

    return mb
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type1::Symbol, input_type2::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_var::Symbol, var_name::Symbol, df::DataFrame;
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            pop_growth::AbstractVector{Float64} = Float64[],
                            shock_name::Nullable{Symbol} = Nullables.Nullable{Symbol}(),
                            density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                            minimize::Bool = false,
                            compute_shockdec_bands::Bool = false)

    # Return only one set of bands if we read in only one draw
    if input_type1 in [:init, :mode, :mean]
        density_bands = [.5]
    end

    # Determine class and product
    class = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast draws
    fcast_series1, transform1 = read_forecast_output(m1, input_type1, cond_type1,
                                                   output_var, var_name, shock_name,
                                                   forecast_string = forecast_string1)
    fcast_series2, transform2 = read_forecast_output(m2, input_type2, cond_type2,
                                                   output_var, var_name, shock_name,
                                                   forecast_string = forecast_string2)

    # Reverse transform
    y0_index = get_y0_index(m1, product)
    yt_index = get_yt_index(m1, product)
    data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,var_name]), NaN))) : fill(NaN, size(df, 1))
    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))
    transformed_series1 = mb_reverse_transform(fcast_series1, transform1, product, class,
                                              y0_index = y0_index, yt_index = yt_index,
                                              data = data,
                                              pop_growth = pop_growth)
    transformed_series2 = mb_reverse_transform(fcast_series2, transform2, product, class,
                                              y0_index = y0_index, yt_index = yt_index,
                                              data = data,
                                              pop_growth = pop_growth)
    @show size(transformed_series1)
    @show size(transformed_series2)
    if size(transformed_series1, 1) < size(transformed_series2, 1)
        mult = Int(size(transformed_series2, 1) / size(transformed_series1, 1))
        transformed_series1 = repeat(transformed_series1, mult, 1)
    else
        mult = Int(size(transformed_series1, 1) / size(transformed_series2, 1))
        transformed_series2 = repeat(transformed_series2, mult, 1)
    end
    @show size(transformed_series1)
    @show size(transformed_series2)

    transformed_series = 0.1*transformed_series1 + 0.9*transformed_series2

    # Compute means and bands
    means = vec(mean(transformed_series, dims= 1))
    bands = if product in [:shockdec, :dettrend, :trend] && !compute_shockdec_bands
        Dict{Symbol,DataFrame}()
    else
        find_density_bands(transformed_series, density_bands, minimize = minimize)
    end
    return means, bands
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_vars::Vector{Symbol};
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            verbose::Symbol = :low, df::DataFrame = DataFrame(),
                            check_empty_columns::Bool = true,
                            kwargs...)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        output_dir = workpath(m1, "forecast")
        println()
        @Base.info "Computing means and bands for input_type = $input_type, cond_type = $cond_type..."
        println("Start time: $(now())")
        println("Means and bands will be saved in $output_dir")
    end
    elapsed_time = @elapsed let
        # Determine full set of output_vars necessary for plotting desired result
        output_vars = add_requisite_output_vars(output_vars)
        if input_type == :prior
            output_vars = setdiff(output_vars, [:bddforecastobs])
        end
        # Load population data and main dataset (required for some transformations)
        if all(var -> get_product(var) == :irf, output_vars)
            population_data, population_forecast = DataFrame(), DataFrame()
        else
            population_data, population_forecast = load_population_growth(m1, verbose = verbose)
            isempty(df) && (df = load_data(m1, check_empty_columns = check_empty_columns, verbose = :none))
        end
        for output_var in output_vars
            prod = get_product(output_var)
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                if prod in [:shockdec, :irf]
                    println("Computing " * string(output_var) * " for shocks:")
                else
                    print("Computing " * string(output_var) * "... ")
                end
            end

            # Compute means and bands
            mb = compute_meansbands(m1, m2, input_type, cond_type1, cond_type2, output_var, df;
                                    forecast_string1 = forecast_string1,
                                    forecast_string2 = forecast_string2,
                                    population_data = population_data,
                                    population_forecast = population_forecast,
                                    verbose = verbose,
                                    kwargs...)
            GC.gc()
        end
    end
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        total_mb_time     = elapsed_time
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                            input_type::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_var::Symbol, df::DataFrame;
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            population_data::DataFrame = DataFrame(),
                            population_forecast::DataFrame = DataFrame(),
                            verbose::Symbol = :none,
                            kwargs...)

    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast metadata
    metadata = get_mb_metadata(m1, input_type, cond_type1, output_var; forecast_string = forecast_string1)

    date_list      = product == :irf ? Date[] : collect(keys(metadata[:date_inds]))
    variable_names = collect(keys(metadata[:indices]))
    pop_growth     = get_mb_population_series(product, population_data, population_forecast, date_list)

    # Compute means and bands
    if product in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                   :bddforecast, :bddforecastut, :bddforecast4q, :dettrend, :trend]
        # Get to work!
        # pmap produces an error for trendobs sometimes, so just doing this iteratively
        mb_vec = Vector{Any}(undef,length(variable_names))
        for i in 1:length(mb_vec)
            mb_vec[i] = compute_meansbands(m1, m2, input_type, cond_type1, cond_type2, output_var,
                                           variable_names[i], df; pop_growth = pop_growth,
                                           forecast_string1 = forecast_string1,
                                           forecast_string2 = forecast_string2,
                                           kwargs...)
        end
        # mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
        #                               pop_growth = pop_growth, forecast_string = forecast_string, kwargs...),
        #               variable_names)

        # Re-assemble pmap outputs
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
            means[!,var_name] = var_means
            bands[var_name] = var_bands
            bands[var_name][!,:date] = date_list
        end

    elseif product in [:shockdec, :irf]
        means = product == :irf ? DataFrame() : DataFrame(date = date_list)
        bands = Dict{Symbol, DataFrame}()

        # Get to work!
        for shock_name in keys(metadata[:shock_indices])
            println(verbose, :high, "  * " * string(shock_name))

            mb_vec = pmap(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
                                          pop_growth = pop_growth, shock_name = Nullables.Nullable(shock_name),
                                          forecast_string1 = forecast_string1,
                                          forecast_string2 = forecast_string2,
                                                         kwargs...),
                          variable_names)

            # Re-assemble pmap outputs
            for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
                means[!, Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_means
                bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_bands
                if product != :irf
                    bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)][!, :date] = date_list
                end
            end
        end

    else
        error("Invalid product: $product")
    end # of if product

    mb = MeansBands(metadata, means, bands)

    # Write to file
    filepath = get_meansbands_output_file(m1, input_type, cond_type1, output_var,
                                          forecast_string = forecast_string1*forecast_string2)
    dirpath = dirname(filepath)
    isdir(dirpath) || mkpath(dirpath)
    JLD2.jldopen(filepath, true, true, true, IOStream) do file
        write(file, "mb", mb)
    end

    sep = prod in [:shockdec, :irf] ? "  " : ""
    println(verbose, :high, sep * "wrote " * basename(filepath))

    return mb
end

function compute_meansbands(m1::AbstractDSGEModel, m2::AbstractDSGEModel, input_type::Symbol,
                            cond_type1::Symbol, cond_type2::Symbol,
                            output_var::Symbol, var_name::Symbol, df::DataFrame;
                            forecast_string1::String = "",
                            forecast_string2::String = "",
                            pop_growth::AbstractVector{Float64} = Float64[],
                            shock_name::Nullable{Symbol} = Nullables.Nullable{Symbol}(),
                            density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                            minimize::Bool = false,
                            compute_shockdec_bands::Bool = false)

    # Return only one set of bands if we read in only one draw
    if input_type in [:init, :mode, :mean]
        density_bands = [.5]
    end

    # Determine class and product
    class = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast draws
    fcast_series1, transform1 = read_forecast_output(m1, input_type, cond_type1,
                                                   output_var, var_name, shock_name,
                                                   forecast_string = forecast_string1)
    fcast_series2, transform2 = read_forecast_output(m2, input_type, cond_type2,
                                                   output_var, var_name, shock_name,
                                                   forecast_string = forecast_string2)

    # Reverse transform
    y0_index = get_y0_index(m1, product)
    yt_index = get_yt_index(m1, product)
    data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,var_name]), NaN))) : fill(NaN, size(df, 1))
    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))
    transformed_series1 = mb_reverse_transform(fcast_series1, transform1, product, class,
                                              y0_index = y0_index, yt_index = yt_index,
                                              data = data,
                                              pop_growth = pop_growth)
    transformed_series2 = mb_reverse_transform(fcast_series2, transform2, product, class,
                                              y0_index = y0_index, yt_index = yt_index,
                                              data = data,
                                              pop_growth = pop_growth)
    transformed_series = 0.5*transformed_series1 + 0.5*transformed_series2


    # Compute means and bands
    means = vec(mean(transformed_series, dims= 1))
    bands = if product in [:shockdec, :dettrend, :trend] && !compute_shockdec_bands
        Dict{Symbol,DataFrame}()
    else
        find_density_bands(transformed_series, density_bands, minimize = minimize)
    end
    return means, bands
end

function mb_reverse_transform(fcast_series::AbstractArray, transform::Function,
                              product::Symbol, class::Symbol;
                              y0_index::Int = -1,
                              data::AbstractVector{Float64} = Float64[],
                              pop_growth::AbstractVector{Float64} = Float64[])
                              # inflation_series::AbstractArray = Float64[])
    # No transformation
    if product in [:histut, :forecastut, :bddforecastut]
        return fcast_series
    end

    if product in [:hist4q, :forecast4q, :bddforecast4q]
        transform4q = get_transform4q(transform)
        use_data = class == :obs && product != :hist4q

        y0s = if use_data && transform4q in [loggrowthtopct_4q_percapita, loggrowthtopct_4q]
            # Sum growth rates y_{t-3}, y_{t-2}, y_{t-1}, and y_t
            data[y0_index+1:y0_index + 3]
        elseif use_data && transform4q in [logleveltopct_4q_percapita, logleveltopct_4q]
            # Divide log levels y_t by y_{t-4}
            data[y0_index:y0_index + 3]
        else
            Float64[]
        end
        reverse_transform(fcast_series, transform4q;
                          fourquarter = true, y0s = y0s,
                          pop_growth = pop_growth)
    # elseif product in [:histlvl, :forecastlvl, :bddforecastlvl]
    #     # NEED TO ADD CASE HERE FOR HISTLVL, ETC# .
        # reverse_transform(fcast_series, transform4q;
        #                   accumulate = true, y0s = y0s,
        #                   pop_growth = pop_growth)
    else
        # Use transformation that doesn't add back population growth for
        # products which are given in deviations
        if product in [:shockdec, :dettrend, :irf]
            transform = get_nopop_transform(transform)
            y0 = 0.0
        else
            y0 = class == :obs ? data[y0_index] : NaN
        end

        reverse_transform(fcast_series, transform;
                          y0 = y0, pop_growth = pop_growth)
    end
end
