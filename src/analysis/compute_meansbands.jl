"""
```
means_bands_all(m, input_type, cond_type, output_vars; forecast_string = "",
    density_bands = [0.5, 0.6, 0.7, 0.8, 0.9], minimize = false, verbose = :low)

means_bands_all(input_type, cond_type, output_vars, input_dir, output_dir,
    filestring_base; forecast_string = "", density_bands = [0.5, 0.6, 0.7, 0.8, 0.9],
    minimize = false, population_mnemonic = Nullable{Symbol}(),
    population_data_file = "", population_forecast_file = "",
    use_population_forecast = true, use_hpfilter = true,
    y0_indexes = Dict{Symbol, Int}(), data = Matrix{Float64}(), verbose = :low)
```

Computes means and bands for pseudo-observables, observables, and shocks, and writes
the results to a file. Two methods are provided. The method that accepts a model
object as an argument uses the model's settings to infer the arguments to the
second method. Users can optionally skip construction of a model object and
manually enter the directory and file names.

### Input Arguments

- `m::AbstractModel`: model object
- `input_type::Symbol`: see `forecast_one`
- `cond_type::Symbol`: see `forecast_one`
- `output_vars::Vector{Symbol}`: see `forecast_one`

**Method 2 only:**

- `input_dir::String`: directory from which the forecast outputs are to
  be read in
- `output_dir::String`: directory to which the means and bands are to be
  saved
- `filestring_base::Vector{String}`: should be equivalent to the result of
  `filestring_base(m)`

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

**Method 2 only:**

- `population_mnemonic::Nullable{Symbol}`: the result of
   `get_setting(m, :population_mnemonic`. Typically `Nullable(:CNP16OV__FRED)`.

- `population_data_file::String`: path to population data (in levels)
  file. In the first method, the following file is used, if it exists:
  `inpath(m, \"raw\", \"population_data_levels_(data_vintage(m)).csv\")`

- `population_forecast_file::S`: path to population forecast (in levels)
  file. In the first method, if `use_population_forecast(m)`, the following file
  is used, if it exists:
  `inpath(m, \"raw\", \"population_forecast_(data_vintage(m)).csv\")`

- `use_population_forecast::Bool`: whether to use population forecast. If
  `false`, use last period of recorded population growth to population adjust
  for forecast periods. Defaults to `true`.

- `use_hpfilter::Bool`: whether to HP filter population data and
  forecast. Defaults to `true`.

- `y0_indexes::Dict{Symbol, Int}`: A `Dict` storing the mapping of products to
  the index of the first period of data needed to compute the product's
  transformation. This is used to compute growth rates and four-quarter
  cumulations (among others). See `get_y0_index`.

- `data::Matrix{Float64}`: pre-loaded `nobs x nperiods` matrix containing the
  (untransformed) data matrix.
"""
function means_bands_all(m::AbstractModel, input_type::Symbol,
                         cond_type::Symbol, output_vars::Vector{Symbol};
                         df::DataFrame = DataFrame(),
                         forecast_string::String = "",
                         density_bands::Vector{Float64} = [0.5, 0.6, 0.7, 0.8, 0.9],
                         minimize::Bool = false,
                         verbose::Symbol = :low)

    ## Step 0: Determine full set of output_vars necessary for plotting desired results
    #          Specifically, if output_vars contains shockdecs but not trend or deterministic trends,
    #          add those

    output_vars = add_requisite_output_vars(output_vars)


    ## Step 1: Population data

    # Parse population mnemonic
    population_mnemonic = parse_population_mnemonic(m)[1]

    # Get filenames for population data and forecast
    if !isnull(population_mnemonic)
        vint = data_vintage(m)
        population_data_file     = inpath(m, "raw", "population_data_levels_$vint.csv")
        population_forecast_file = inpath(m, "raw", "population_forecast_$vint.csv")
    else
        population_data_file = ""
        population_forecast_file = ""
    end

    ## Step 2: Load main dataset (required for some transformations),
    ##         specify which period to use as first t-1 period for computing
    ##         growth rates, and which period is the last historical period

    # Load dataset if not provided
    if isempty(df)
        df = load_data(m, verbose = :none)
    end
    data = df_to_matrix(m, df)

    # Specify the t-1 period for each product
    y0_indexes = Dict{Symbol,Int}()
    products   = unique(map(get_product, output_vars))
    for prod in products
        y0_indexes[prod] = get_y0_index(m, prod)
    end

    ## Step 3: Fetch the model string and the input and output directories we'll read from/write to

    input_dir  = rawpath(m, "forecast")
    output_dir = workpath(m, "forecast")
    base       = filestring_base(m)

    ## Step 4: We have everything we need; appeal to model-object-agnostic function

    means_bands_all(input_type, cond_type, output_vars,
                    input_dir, output_dir, base;
                    forecast_string = forecast_string,
                    density_bands = density_bands, minimize = minimize,
                    population_mnemonic = population_mnemonic,
                    population_data_file = population_data_file,
                    population_forecast_file = population_forecast_file,
                    use_population_forecast = use_population_forecast(m),
                    use_hpfilter = hpfilter_population(m),
                    y0_indexes = y0_indexes, data = data,
                    compute_shockdec_bands = get_setting(m, :compute_shockdec_bands),
                    verbose = verbose)
end

function means_bands_all(input_type::Symbol, cond_type::Symbol, output_vars::Vector{Symbol},
                         input_dir::String, output_dir::String,
                         filestring_base::Vector{String};
                         forecast_string::String = "",
                         density_bands::Vector{Float64} = [0.5, 0.6, 0.7, 0.8, 0.9],
                         minimize::Bool = false,
                         population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                         population_data_file::String = "",
                         population_forecast_file::String = "",
                         use_population_forecast::Bool = true,
                         use_hpfilter::Bool = true,
                         y0_indexes::Dict{Symbol,Int} = Dict{Symbol,Int}(),
                         data = Matrix{Float64}(),
                         compute_shockdec_bands::Bool = false,
                         verbose::Symbol = :low)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println()
        info("Computing means and bands for input_type = $input_type, cond_type = $cond_type...")
        println("Start time: $(now())")
        println("Means and bands will be saved in $output_dir")
        tic()
    end

    ## Step 1: Filter population history and forecast and compute growth rates

    dlfiltered_population_data, dlfiltered_population_forecast =
        if isnull(population_mnemonic)
            DataFrame(), DataFrame()
        else
            load_population_growth(population_data_file, population_forecast_file,
                                   get(population_mnemonic);
                                   use_population_forecast = use_population_forecast,
                                   use_hpfilter = use_hpfilter,
                                   verbose = verbose)
        end

    ## Step 2: Set up filenames for MeansBands input and output files.
    # MeansBands output filenames are the same as forecast output filenames, but with an "mb" prefix.

    mb_input_files = get_meansbands_input_files(input_dir, filestring_base,
                                                input_type, cond_type, output_vars;
                                                forecast_string = forecast_string)
    mb_output_files = get_meansbands_output_files(output_dir, filestring_base,
                                                  input_type, cond_type, output_vars,
                                                  forecast_string = forecast_string)

    ## Step 3: Compute means and bands for each output variable, and write to a file.

    for output_var in output_vars

        # Which product are we dealing with? Need this to index out of y0_indexes.
        prod     = get_product(output_var)
        y0_index = isempty(y0_indexes[prod]) ? -1 : y0_indexes[prod]

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            if prod in [:shockdec, :irf]
                println("Computing " * string(output_var) * " for shocks:")
            else
                print("Computing " * string(output_var) * "... ")
            end
        end

        # compute means and bands object
        mb = means_bands(input_type, cond_type, output_var,
                         mb_input_files, density_bands = density_bands, minimize = minimize,
                         forecast_string = forecast_string,
                         population_data = dlfiltered_population_data,
                         population_forecast = dlfiltered_population_forecast,
                         population_mnemonic = Nullable(:population_growth),
                         y0_index = y0_index, data = data, verbose = verbose,
                         compute_shockdec_bands = compute_shockdec_bands)

        # write to file
        filepath = mb_output_files[output_var]
        dirpath  = dirname(filepath)
        isdir(dirpath) || mkpath(dirpath)
        jldopen(filepath, "w") do file
            write(file, "mb", mb)
        end
        gc()

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            sep = prod in [:shockdec, :irf] ? "  " : ""
            println(sep * "wrote " * basename(filepath))
        end
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        total_mb_time     = toq()
        total_mb_time_min = total_mb_time/60

        println("\nTotal time to compute means and bands: " * string(total_mb_time_min) * " minutes")
        println("Computation of means and bands complete: " * string(now()))
    end
end

"""
```
means_bands(input_type, cond_type, output_var,
    meansbands_input_files::Dict{Symbol, String}; forecast_string = "",
    density_bands = [0.5, 0.6, 0.7, 0.8, 0.9], minimize = false,
    population_mnemonic = Nullable{Symbol}(), population_data = DataFrame(),
    population_forecast = DataFrame(), y0_index = -1, data = Matrix{Float64}(),
    verbose = :low)
```

Computes means and bands for a single `output_var`.
"""
function means_bands(input_type::Symbol,
                     cond_type::Symbol,
                     output_var::Symbol,
                     meansbands_input_files::Dict{Symbol, String};
                     forecast_string::String = "",
                     density_bands::Vector{Float64} = [0.5, 0.6, 0.7, 0.8, 0.9],
                     minimize::Bool = false,
                     population_data::DataFrame = DataFrame(),
                     population_forecast::DataFrame = DataFrame(),
                     population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                     y0_index::Int = -1,
                     data::Matrix{Float64} = Matrix{Float64}(),
                     verbose::Symbol = :none,
                     compute_shockdec_bands::Bool = false)

    # Return only one set of bands if we read in only one draw
    if input_type in [:init, :mode, :mean]
        density_bands = [.5]
    end

    ## Step 1: Determine the class of variable we are working with (pseudos?
    ##         observables? etc) and the product we are computing (forecast?
    ##         history? shockdec?)

    class   = get_class(output_var)
    product = get_product(output_var)

    ## Step 2: Read in raw forecast output and metadata (transformations,
    ##         mappings from Symbols to indices, etc)

    forecast_output_file = meansbands_input_files[output_var]
    metadata, mb_metadata =
        get_mb_metadata(input_type, cond_type, output_var, forecast_output_file;
                        forecast_string = forecast_string)

    date_list         = product == :irf ? Vector{Date}() : collect(keys(mb_metadata[:date_inds]))
    variable_names    = collect(keys(mb_metadata[:indices]))
    population_series = if isnull(population_mnemonic)
        Vector{Float64}()
    else
        get_mb_population_series(product, get(population_mnemonic), population_data,
                                 population_forecast, date_list)
    end

    ## Step 3: Compute means and bands

    if product in [:hist, :forecast, :dettrend, :trend, :hist4q, :forecast4q, :bddforecast, :bddforecast4q]

        # Get to work!
        mb_vec = pmap(var_name -> compute_means_bands(class, product, var_name, forecast_output_file;
                          data = data, population_series = population_series, y0_index = y0_index,
                          density_bands = density_bands, minimize = minimize,
                          compute_shockdec_bands = compute_shockdec_bands),
                      variable_names)

        # Re-assemble pmap outputs
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()
        for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
            means[var_name] = var_means
            bands[var_name] = var_bands
            bands[var_name][:date] = date_list
        end

    elseif product in [:shockdec, :irf]

        means = product == :irf ? DataFrame() : DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        mb_metadata[:shock_indices] = metadata[:shock_indices]

        # Get to work!
        for shock_name in keys(mb_metadata[:shock_indices])
            if VERBOSITY[verbose] >= VERBOSITY[:high]
                println("  * " * string(shock_name))
            end

            mb_vec = pmap(var_name -> compute_means_bands(class, product, var_name, forecast_output_file;
                              data = data, population_series = population_series, y0_index = y0_index,
                              shock_name = Nullable(shock_name), density_bands = density_bands,
                              minimize = minimize, compute_shockdec_bands = compute_shockdec_bands),
                          variable_names)

            # Re-assemble pmap outputs
            for (var_name, (var_means, var_bands)) in zip(variable_names, mb_vec)
                means[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_means
                bands[Symbol(var_name, DSGE_SHOCKDEC_DELIM, shock_name)] = var_bands
                if product != :irf
                    bands[Symbol(var_name,DSGE_SHOCKDEC_DELIM, shock_name)][:date] = date_list
                end
            end
        end
    end # of if product

    return MeansBands(mb_metadata, means, bands)
end

function compute_means_bands(class::Symbol,
                             product::Symbol,
                             var_name::Symbol,
                             filename::String;
                             data::Matrix{Float64} = Matrix{Float64}(),
                             population_series::Vector{Float64} = Vector{Float64}(),
                             y0_index::Int = -1,
                             shock_name::Nullable{Symbol} = Nullable{Symbol}(),
                             density_bands::Vector{Float64} = [0.5,0.6,0.7,0.8,0.9],
                             minimize::Bool = false,
                             compute_shockdec_bands::Bool = false)

    fcast_series, transform, var_ind, date_list = jldopen(filename, "r") do file
        # Read forecast output
        fcast_series = if isnull(shock_name)
            read_forecast_output(file, class, product, var_name)
        else
            read_forecast_output(file, class, product, var_name, get(shock_name))
        end

        # Parse transform
        class_long = get_class_longname(class)
        transforms = read(file, string(class_long) * "_revtransforms")
        transform = parse_transform(transforms[var_name])

        # Get variable index
        indices = read(file, string(class_long) * "_indices")
        var_ind = indices[var_name]

        # Read date list
        date_list = if product != :irf
            collect(keys(read(file, "date_indices")))
        else
            Date[]
        end

        fcast_series, transform, var_ind, date_list
    end

    # The `fcast_output` for trends only is of size `ndraws` x `nvars`. We
    # need to use `repeat` below because population adjustments will be
    # different in each period. Now we have something of size `ndraws` x
    # `nvars` x `nperiods`
    if product == :trend
        fcast_series = repeat(fcast_series, outer = [1, length(date_list)])
    end

    # Do we want to use the data for y0 and y0s?
    use_data = class == :obs && !(product in [:irf, :hist4q])

    # Reverse transform
    if product in [:hist4q, :forecast4q, :bddforecast4q]
        transform4q = get_transform4q(transform)

        y0s = if transform4q in [loggrowthtopct_4q_percapita, loggrowthtopct_4q]
            # Sum growth rates y_{t-3}, y_{t-2}, y_{t-1}, and y_t
            use_data ? data[var_ind, y0_index+1:end] : fill(NaN, 3)
        elseif transform4q in [logleveltopct_4q_percapita, logleveltopct_4q]
            # Divide log levels y_t by y_{t-4}
            use_data ? data[var_ind, y0_index:end] : fill(NaN, 4)
        else
            Float64[]
        end

        transformed_series = reverse_transform(fcast_series, transform4q;
                                               fourquarter = true,
                                               y0s = y0s, pop_growth = population_series)
    else
        # Use IRF transform if necessary
        if product == :irf
            transform = get_irf_transform(transform)
        end

        y0 = use_data ? data[var_ind, y0_index] : NaN
        transformed_series = reverse_transform(fcast_series, transform;
                                               y0 = y0, pop_growth = population_series)
    end

    # Compute means and bands of transformed series
    means = vec(mean(transformed_series, 1))
    bands = if product in [:shockdec, :dettrend, :trend] && compute_shockdec_bands
        Dict{Symbol,DataFrame}()
    else
        find_density_bands(transformed_series, density_bands, minimize = minimize)
    end
    return means, bands
end
