"""
```
meansbands_arithmetic(op, m1, m2, input_type, cond_type, output_vars; forecast_string = "",
    verbose = :low, kwargs...)

meansbands_arithmetic(op, m1, m2, input_type, cond_type, output_var, df; forecast_string = "",
    population_data = DataFrame(), population_forecast = DataFrame(),
    verbose = :none, kwargs...)

meansbands_arithmetic(op, m1, m2, input_type, cond_type, output_var, var_name, df;
    forecast_string = "", population_data = DataFrame(),
    population_forecast = DataFrame(), verbose = :low, kwargs...)
```

Calculates arithmetic on means and bands for pseudo-observables, observables, and shocks
from two different computations, and write
the results to a file. Other methods are for one `output_var` and one `var_name` respectively.
For additional keywords, see the list below.

### Inputs
- `op::Function`: arithmetic operation. Must be +, -, *, /, or ^.

- `m1::AbstractDSGEModel`: DSGE model object corresponding to the first `MeansBands`

- `m2::AbstractDSGEModel`:DSGE model object corresponding to the second `MeansBands`

- `cond_type::Symbol`: type of conditional forecast.

- `input_type::Symbol`: type of forecast used. Typically `:full` or `:mode`.

### Keyword Arguments

- `forecast_string1::String`: forecast identifier string (the value
  \"fcid=value\" in the forecast output filename) for `m1`.

- `forecast_string2::String`: forecast identifier string (the value
  \"fcid=value\" in the forecast output filename) for `m2`.

- `density_bands::Vector{Float64}`: a vector of percent values (between 0 and 1) for
  which to compute density bands

- `df1::DataFrame1`: data to be used when computing the `MeansBands` for model `m1`.

- `df2::DataFrame2`: data to be used when computing the `MeansBands` for model `m2`.

- `minimize::Bool`: if `true`, choose shortest interval, otherwise just chop off
  lowest and highest (percent/2)

- `verbose`: level of error messages to be printed to screen. One of `:none`,
  `:low`, `:high`

- `check_empty_columns::Bool`: if true, throw an error if
    calling load_data yields an empty column.

- `do_cond_obs_shocks1::Bool`: if true, append "cond_obs_shocks"
    to names of computed `MeansBands` files corresponding
    to model `m1`.

- `do_cond_obs_shocks2::Bool`: if true, append "cond_obs_shocks"
    to names of computed `MeansBands` files corresponding
    to model `m2`.

- `output_model::Int`: specifies to which model's saveroot
    output will be written. Must be either 1 or 2.
```
"""
function meansbands_arithmetic(op::Function, m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                               input_type::Symbol,
                               cond_type::Symbol, output_var::Symbol;
                               forecast_string1::String = "",
                               forecast_string2::String = "",
                               verbose::Symbol = :low, df1::DataFrame = DataFrame(),
                               df2::DataFrame = DataFrame(),
                               check_empty_columns::Bool = true,
                               do_cond_obs_shocks1::Bool = false,
                               do_cond_obs_shocks2::Bool = false,
                               output_model::Int = 1,
                               kwargs...)
    mb = meansbands_arithmetic(op, m1, m2, input_type, cond_type, [output_var];
                               forecast_string1 = forecast_string1,
                               forecast_string2 = forecast_string2,
                               verbose = verbose, df1 = df1,
                               df2 = df2, check_empty_columns = check_empty_columns,
                               do_cond_obs_shocks1 = do_cond_obs_shocks1,
                               do_cond_obs_shocks2 = do_cond_obs_shocks2,
                               output_model = output_model,
                               kwargs...)
    return mb
end
function meansbands_arithmetic(op::Function, m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                               input_type::Symbol,
                               cond_type::Symbol, output_vars::Vector{Symbol};
                               forecast_string1::String = "",
                               forecast_string2::String = "",
                               verbose::Symbol = :low, df1::DataFrame = DataFrame(),
                               df2::DataFrame = DataFrame(),
                               check_empty_columns::Bool = true,
                               do_cond_obs_shocks1::Bool = false,
                               do_cond_obs_shocks2::Bool = false,
                               output_model::Int = 1,
                               kwargs...)
    if !(op in [Base.:+, Base.:-, Base.:*, Base.:/, Base.:^])
        error("Operator $op is not an accepted arithmetic operator. The only allowed operators are" *
              " +, -, *, /, and ^.")
    end
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        output_dir = workpath(m1, "forecast")
        println()
        @Base.info "Computing arithmetic on means and bands for input_type = $input_type, cond_type = $cond_type..."
        println("Start time: $(now())")
        # println("Means and bands will be saved in $output_dir")
    end

    # Calculate MeansBands
    mb_dict = Dict{Symbol,MeansBands}()
    elapsed_time = @elapsed let
        # Determine full set of output_vars necessary for plotting desired result
        output_vars = add_requisite_output_vars(output_vars)
        if input_type == :prior
            output_vars = setdiff(output_vars, [:bddforecastobs])
        end
        # Load population data and main dataset (required for some transformations)
        if all(var -> get_product(var) == :irf, output_vars)
            population_data1, population_forecast1 = DataFrame(), DataFrame()
            population_data2, population_forecast2 = DataFrame(), DataFrame()
        else
            population_data1, population_forecast1 = load_population_growth(m1, verbose = verbose)
            isempty(df1) && (df1 = load_data(m1, check_empty_columns = check_empty_columns, verbose = :none))
            population_data2, population_forecast2 = load_population_growth(m2, verbose = verbose)
            isempty(df2) && (df2 = load_data(m2, check_empty_columns = check_empty_columns, verbose = :none))
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
            mb_dict[output_var] = meansbands_arithmetic(op, m1, m2, input_type,
                                                        cond_type, output_var, df1, df2;
                                                        forecast_string1 = forecast_string1,
                                                        forecast_string2 = forecast_string2,
                                                        population_data1 = population_data1,
                                                        population_forecast1 = population_forecast1,
                                                        population_data2 = population_data2,
                                                        population_forecast2 = population_forecast2,
                                                        verbose = verbose,
                                                        do_cond_obs_shocks1 = do_cond_obs_shocks1,
                                                        do_cond_obs_shocks2 = do_cond_obs_shocks2,
                                                        output_model = output_model,
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
    return mb_dict
end

function meansbands_arithmetic(op::Function, m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                               input_type::Symbol, cond_type::Symbol,
                               output_var::Symbol, df1::DataFrame, df2::DataFrame;
                               forecast_string1::String = "",
                               forecast_string2::String = "",
                               population_data1::DataFrame = DataFrame(),
                               population_forecast1::DataFrame = DataFrame(),
                               population_data2::DataFrame = DataFrame(),
                               population_forecast2::DataFrame = DataFrame(),
                               verbose::Symbol = :none,
                               do_cond_obs_shocks1::Bool = false,
                               do_cond_obs_shocks2::Bool = false,
                               output_model::Int = 1,
                               kwargs...)
    if output_model != 1 && output_model != 2
        error("The keyword `output_model` can be only 1 or 2")
    end

    # Determine class and product
    class   = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast metadata
    metadata1 = get_mb_metadata(m1, input_type, cond_type, output_var; forecast_string = forecast_string1)
    metadata2 = get_mb_metadata(m2, input_type, cond_type, output_var; forecast_string = forecast_string2)

    date_list1      = product == :irf ? Date[] : collect(keys(metadata1[:date_inds]))
    date_list2      = product == :irf ? Date[] : collect(keys(metadata2[:date_inds]))
    variable_names1 = collect(keys(metadata1[:indices]))
    variable_names2 = collect(keys(metadata2[:indices]))
    if !isempty(setdiff(date_list1, date_list2))
        error("Date indices of forecasts for m1 and m2 for input_type $input_type, " *
              "cond_type $cond_type, and output_var $output_var do not match")
    end
    if !isempty(setdiff(variable_names1, variable_names2))
        error("Variable names from forecasts for m1 and m2 for input_type $input_type, " *
              "cond_type $cond_type, and output_var $output_var do not match")
    end
    pop_growth1     = get_mb_population_series(product, population_data1, population_forecast1, date_list1)
    pop_growth2     = get_mb_population_series(product, population_data2, population_forecast2, date_list2)

    # Compute means and bands
    if product in [:hist, :histut, :hist4q, :forecast, :forecastut, :forecast4q,
                   :bddforecast, :bddforecastut, :bddforecast4q, :dettrend, :trend]
        # Get to work!
        # pmap produces an error for trendobs sometimes, so just doing this iteratively
        mb_vec = Vector{Any}(undef,length(variable_names1))
        for i in 1:length(mb_vec)
            mb_vec[i] = meansbands_arithmetic(op, m1, m2, input_type, cond_type, output_var,
                                              variable_names1[i], df1, df2; pop_growth1 = pop_growth1,
                                              pop_growth2 = pop_growth2,
                                              forecast_string1 = forecast_string1,
                                              forecast_string2 = forecast_string2,
                                              do_cond_obs_shocks1 = do_cond_obs_shocks1,
                                              do_cond_obs_shocks2 = do_cond_obs_shocks2,
                                              kwargs...)
        end

        # Re-assemble pmap outputs
        means = DataFrame(date = date_list1)
        bands = Dict{Symbol,DataFrame}()

        for (var_name, (var_means, var_bands)) in zip(variable_names1, mb_vec)
            means[!,var_name] = var_means
            bands[var_name] = var_bands
            bands[var_name][!,:date] = date_list1
        end

    elseif product in [:shockdec, :irf]
        means = product == :irf ? DataFrame() : DataFrame(date = date_list1)
        bands = Dict{Symbol, DataFrame}()

        # Get to work!
        for shock_name in keys(metadata[:shock_indices])
            println(verbose, :high, "  * " * string(shock_name))

            mb_vec = pmap(var_name -> meansbands_arithmetic(op, m1, m2, input_type, cond_type, output_var,
                                                            var_name, df1, df2;
                                                            pop_growth1 = pop_growth1,
                                                            pop_growth2 = pop_growth2,
                                                            shock_name = Nullables.Nullable(shock_name),
                                                            forecast_string1 = forecast_string1,
                                                            forecast_string2 = forecast_string2,
                                                            do_cond_obs_shocks1 = do_cond_obs_shocks1,
                                                            do_cond_obs_shocks2 = do_cond_obs_shocks2,
                                                            kwargs...),
                          variable_names1)

            # Re-assemble pmap outputs
            for (var_name, (var_means, var_bands)) in zip(variable_names1, mb_vec)
                means[!, Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_means
                bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_bands
                if product != :irf
                    bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)][!, :date] = date_list1
                end
            end
        end

    else
        error("Invalid product: $product")
    end # of if product

    mb = MeansBands(metadata1, means, bands)
    filepath =
        get_meansbands_arithmetic_output_file(op, m1, m2,
                                              input_type, cond_type,
                                              output_var,
                                              output_model = output_model,
                                              forecast_string1 =
                                              forecast_string1,
                                              forecast_string2 =
                                              forecast_string2,
                                              do_cond_obs_shocks1 =
                                              do_cond_obs_shocks1,
                                              do_cond_obs_shocks2 =
                                              do_cond_obs_shocks2)

    dirpath = dirname(filepath)
    isdir(dirpath) || mkpath(dirpath)
    JLD2.jldopen(filepath, true, true, true, IOStream) do file
        write(file, "mb", mb)
    end

    sep = prod in [:shockdec, :irf] ? "  " : ""
    println(verbose, :high, sep * "wrote " * basename(filepath))

    return mb
end

function meansbands_arithmetic(op::Function, m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                               input_type::Symbol, cond_type::Symbol,
                               output_var::Symbol, var_name::Symbol,
                               df1::DataFrame, df2::DataFrame;
                               forecast_string1::String = "",
                               forecast_string2::String = "",
                               pop_growth1::AbstractVector{Float64} = Float64[],
                               pop_growth2::AbstractVector{Float64} = Float64[],
                               shock_name::Nullable{Symbol} = Nullables.Nullable{Symbol}(),
                               density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                               minimize::Bool = false,
                               compute_shockdec_bands::Bool = false,
                               do_cond_obs_shocks1::Bool = false,
                               do_cond_obs_shocks2::Bool = false)

    # Return only one set of bands if we read in only one draw
    if input_type in [:init, :mode, :mean]
        density_bands = [.5]
    end

    # Determine class and product
    class = get_class(output_var)
    product = get_product(output_var)

    # Read in forecast draws
    fcast_series1, transform1 = read_forecast_output(m1, input_type, cond_type,
                                                     output_var, var_name, shock_name,
                                                     forecast_string = forecast_string1,
                                                     do_cond_obs_shocks = do_cond_obs_shocks1)
    fcast_series2, transform2 = read_forecast_output(m2, input_type, cond_type,
                                                     output_var, var_name, shock_name,
                                                     forecast_string = forecast_string2,
                                                     do_cond_obs_shocks = do_cond_obs_shocks2)

    # Reverse transform
    y0_index1 = get_y0_index(m1, product)
    y0_index2 = get_y0_index(m2, product)
    data1 = class == :obs && product != :irf ? Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df1[!,var_name]), NaN))) : fill(NaN, size(df1, 1))
    data2 = class == :obs && product != :irf ? Float64.(collect(Missings.replace(Vector{Union{Missing, Float64}}(df2[!,var_name]), NaN))) : fill(NaN, size(df2, 1))
    # data = class == :obs && product != :irf ? Float64.(collect(Missings.replace(df[:,var_name], NaN))) : fill(NaN, size(df, 1))
    transformed_series1 = mb_reverse_transform(fcast_series1, transform1, product, class,
                                               y0_index = y0_index1, data = data1,
                                               pop_growth = pop_growth1)
    transformed_series2 = mb_reverse_transform(fcast_series2, transform2, product, class,
                                               y0_index = y0_index2, data = data2,
                                               pop_growth = pop_growth2)

    # Compute means and bands
    if op in [Base.:+, Base.:-]
        transformed_series = op(transformed_series1, transformed_series2)
    elseif op in [Base.:*, Base.:/, Base.:^]
        transformed_series = map((x,y) -> op(x,y), transformed_series1, transformed_series2)
    else
        error("Operator $op is not an arithmetic operator. It must be one of " *
              "+, -, /, *, or ^")
    end
    means = vec(mean(transformed_series, dims= 1))
    bands = if product in [:shockdec, :dettrend, :trend] && !compute_shockdec_bands
        Dict{Symbol,DataFrame}()
    else
        find_density_bands(transformed_series, density_bands, minimize = minimize)
    end
    return means, bands
end
