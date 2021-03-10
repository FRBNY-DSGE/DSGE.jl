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

- `forecast_string::String = ""`: forecast identifier string (the value
  \"fcid=value\" in the forecast output filename). Required when
  `input_type == :subset`

- `density_bands::Vector{Float64} = [.5, .6, .7, .8 .9]`: a vector of percent values (between 0 and 1) for
  which to compute density bands.

- `minimize::Bool = false`: if `true`, choose shortest interval, otherwise just chop off
  lowest and highest (percent/2)

- `verbose::Symbol = :low`: level of error messages to be printed to screen. One of `:none`,
  `:low`, `:high`

- `bdd_fcast::Bool = true`: if true, calculate bounded forecasts

- `skipnan::Bool = false`: if true, remove any NaNs found in the raw forecast output series

- `df::DataFrame = DataFrame()`: if an empty DataFrame, then the function will attempt to
    load the data using `load_data`.

- `check_empty_columns::Bool = true`: if true, throw an error if
    calling load_data yields an empty column.

- `pseudo2data::AbstractDict{Symbol, Symbol} = Dict()`: Maps the names of pseudo-observables
    which require historical data during transformation to the name of the required historical series.
"""
function compute_meansbands(m::AbstractDSGEModel, input_type::Symbol,
                            cond_type::Symbol, output_vars::Vector{Symbol};
                            forecast_string::String = "",
                            verbose::Symbol = :low, df::DataFrame = DataFrame(),
                            check_empty_columns::Bool = true,
                            bdd_fcast::Bool = true, skipnan::Bool = false,
                            pseudo2data::AbstractDict{Symbol, Symbol} = Dict{Symbol, Symbol}(),
                            variable_names::Vector{Symbol} = Vector{Symbol}(undef, 0),
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
                                    skipnan = skipnan, pseudo2data = pseudo2data,
                                    variable_names = variable_names,
                                    verbose = verbose,
                                    bdd_fcast = bdd_fcast,
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
                            skipnan::Bool = false,
                            pseudo2data::AbstractDict{Symbol, Symbol} = Dict{Symbol, Symbol}(),
                            variable_names::Vector{Symbol} = Vector{Symbol}(undef, 0),
                            verbose::Symbol = :none,
                            bdd_fcast::Bool = true,
                            kwargs...)

    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast metadata
    metadata = get_mb_metadata(m, input_type, cond_type, output_var; forecast_string = forecast_string)

    date_list      = product == :irf ? Date[] : collect(keys(metadata[:date_inds]))
    if isempty(variable_names)
        variable_names = collect(keys(metadata[:indices]))
    end
    pop_growth     = get_mb_population_series(product, population_data, population_forecast, date_list)

    # Compute means and bands
    map_fcn = get_setting(m, :use_parallel_workers) ? pmap : map
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
                                               skipnan = skipnan, pseudo2data = pseudo2data,
                                               bdd_fcast = bdd_fcast,
                                               kwargs...)
            end
        else
            mb_vec = map_fcn(var_name -> compute_meansbands(m, input_type, cond_type, output_var, var_name, df;
                                                            pop_growth = pop_growth, forecast_string = forecast_string,
                                                            pseudo2data = pseudo2data,
                                                            bdd_fcast = bdd_fcast,
                                                            skipnan = skipnan, kwargs...),
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
                                                         skipnan = skipnan, pseudo2data = pseudo2data, bdd_fcast = bdd_fcast, kwargs...),
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
                            skipnan::Bool = false,
                            shock_name::Nullable{Symbol} = Nullables.Nullable{Symbol}(),
                            density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                            minimize::Bool = false,
                            pseudo2data::AbstractDict{Symbol, Symbol} = Dict{Symbol, Symbol}(),
                            compute_shockdec_bands::Bool = false,
                            bdd_fcast::Bool = true)

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
    data = if haskey(pseudo2data, var_name) && product != :irf
        Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,pseudo2data[var_name]]), NaN)))
    elseif class == :obs && product != :irf
        Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,var_name]), NaN)))
    else
        fill(NaN, size(df, 1))
    end

    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))
    transformed_series = mb_reverse_transform(fcast_series, transform, product, class,
                                              y0_index = y0_index, data = data,
                                              pop_growth = pop_growth, use_data = haskey(pseudo2data, var_name))

    replace!(transformed_series, missing=>NaN)

    # Handle NaNs
    if skipnan && !(output_var in [:histobs, :hist4qobs]) && any(isnan.(transformed_series))
        # Remove rows with NaNs
        nanrows = vec(mapslices(x -> all(isnan.(x)), transformed_series, dims = Int[2]))
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

"""
```
compute_meansbands(models, input_types, cond_types, output_vars; forecast_strings = [],
    verbose = :low, kwargs...)

compute_meansbands(models, input_types, cond_types, output_vars, df; forecast_strings = [],
    population_data = DataFrame(), population_forecast = DataFrame(),
    verbose = :none, kwargs...)

compute_meansbands(models, input_types, cond_types, output_var, var_name, df;
    forecast_strings = [], population_data = DataFrame(),
    population_forecast = DataFrame(), verbose = :low,
    kwargs...)
```

Compute means and bands for pseudo-observables, observables, and shocks from multiple forecasts, and write
the results to a file. Other methods are for one `output_var` and one `var_name` respectively.

### Keyword Arguments

- `forecast_strings::Vector{String} = []`: forecast identifier strings for each forecast (the value
  \"fcid=value\" in the forecast output filename for individual forecasts). Required when
  `input_type == :subset`

- `density_bands::Vector{Float64} = [.5, .6, .7, .8 .9]`: a vector of percent values (between 0 and 1) for
  which to compute density bands.

- `minimize::Bool = false`: if `true`, choose shortest interval, otherwise just chop off
  lowest and highest (percent/2)

- `verbose::Symbol = :low`: level of error messages to be printed to screen. One of `:none`,
  `:low`, `:high`

- `bdd_fcast::Bool = true`: if true, calculate bounded forecasts

- `skipnan::Bool = false`: if true, remove any NaNs found in the raw forecast output series

- `ndraws::Int = 20000`: number of draws from the underlying scenarios

- `df::DataFrame = DataFrame()`: if an empty DataFrame, then the function will attempt to
    load the data using `load_data`.

- `check_empty_columns::Bool = true`: if true, throw an error if
    calling load_data yields an empty column.

- `pseudo2data::AbstractDict{Symbol, Symbol} = Dict()`: Maps the names of pseudo-observables
    which require historical data during transformation to the name of the required historical series.
"""
function compute_meansbands(models::Vector{<: AbstractDSGEModel},
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
                            skipnan::Bool = false, ndraws::Int = 20000,
                            pseudo2data::AbstractDict{Symbol, Symbol} = Dict{Symbol, Symbol}(),
                            verbose::Symbol = :low,
                            bdd_fcast::Bool = true,
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
        output_vars = add_requisite_output_vars(output_vars, bdd_fcast=bdd_fcast)
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
                                    skipnan = skipnan, ndraws = ndraws,
                                    pseudo2data = pseudo2data,
                                    verbose = verbose,
                                    bdd_fcast = bdd_fcast,
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

function compute_meansbands(models::Vector{<: AbstractDSGEModel},
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
                            skipnan::Bool = false, ndraws::Int = 20000,
                            pseudo2data::AbstractDict{Symbol, Symbol} = Dict{Symbol, Symbol}(),
                            verbose::Symbol = :none,
                            bdd_fcast::Bool = true,
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
                                               weights = weights, skipnan = skipnan,
                                               pseudo2data = pseudo2data,
                                               ndraws = ndraws, bdd_fcast = bdd_fcast,
                                               kwargs...)
            end
        else
            mb_vec = pmap(var_name -> compute_meansbands(models, input_types, cond_types, output_var, var_name, df;
                                                         pop_growth = pop_growth,
                                                         forecast_strings = forecast_strings, weights = weights,
                                                         skipnan = skipnan, ndraws = ndraws,
                                                         pseudo2data = pseudo2data,
                                                         bdd_fcast = bdd_fcast, kwargs...),
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
                                                         skipnan = skipnan, ndraws = ndraws,
                                                         pseudo2data = pseudo2data, kwargs...),
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
                            minimize::Bool = false, skipnan::Bool = false,
                            compute_shockdec_bands::Bool = false, ndraws::Int = 20000,
                            pseudo2data::AbstractDict{Symbol, Symbol} = Dict{Symbol, Symbol}(),
                            bdd_fcast::Bool = true)

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
    data = if haskey(pseudo2data, var_name) && product != :irf
        Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,pseudo2data[var_name]]), NaN)))
    elseif class == :obs && product != :irf
        Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df[!,var_name]), NaN)))
    else
        fill(NaN, size(df, 1))
    end
    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))


    fcast_series = Matrix{Float64}(undef, 0, size(fcast_seriess[1], 2))
    if all(map(x->size(x, 1) > 1, fcast_seriess))
        for i in 1:length(models)
            sampled = sample(1:size(fcast_seriess[i], 1), Int(ndraws * weights[i]))
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
                                              pop_growth = pop_growth,  use_data = haskey(pseudo2data, var_name))

    replace!(transformed_series, missing=>NaN)

    # Handle NaNs
    if skipnan && !(output_var in [:histobs, :hist4qobs]) && (ismissing(any(isnan.(transformed_series))) ? true : any(isnan.(transformed_series)))
    # if skipnan && !(output_var in [:histobs, :hist4qobs]) && any(isnan.(transformed_series)) # TODO: figure out why missings occur
        # Remove rows with NaNs
        nanrows = vec(mapslices(x -> ismissing(all(isnan.(x))) ? all(ismissing.(x)) : all(isnan.(x)), transformed_series, dims = Int[2]))
        transformed_series = transformed_series[.!nanrows, :]
    end

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
                              pop_growth::AbstractVector{Float64} = Float64[],
                              use_data::Bool = false)
                              # inflation_series::AbstractArray = Float64[])
    # No transformation
    if product in [:histut, :forecastut, :bddforecastut]
        return fcast_series
    end

    if product in [:hist4q, :forecast4q, :bddforecast4q]
        transform4q = get_transform4q(transform)
        use_data = use_data ? true : (class == :obs && product != :hist4q)

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
