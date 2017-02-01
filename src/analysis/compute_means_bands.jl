"""
```
means_bands_all(m, input_type, cond_type, output_vars;
                density_bands = [0.5, 0.6, 0.7, 0.8, 0.9], minimize = false,
                forecast_string = "", load_dataset = true, load_population_data = true,
                population_forecast_file = "", verbose :low)


means_bands_all(input_type, cond_type, output_vars, meansbands_input_files;
                density_bands = [0.5, 0.6, 0.7, 0.8, 0.9], minimize = false,
                forecast_string = "", output_dir = "", population_data = DataFrame(),
                population_mnemonic = Nullable{Symbol}(),
                population_forecast_file = "", y0_indexes = Dict{Symbol,Nullable{Int}}(),
                data = Matrix{T}(), verbose::Symbol = :low)
```

Computes means and bands for pseudoobservables and observables, and writes
results to a file. Two methods are provided. The method that accepts a model
object as an argument uses the model's settings to infer `forecast_files`; then
appeals to the second method. Users can optionally skip construction of a model
object and manually enter `forecast_files`.

Below, `T<:AbstractFloat` and `S<:AbstractString`:

### Input Arguments

- `m`: model object
- `input_type::Symbol`: see `forecast_all`
- `cond_type::Symbol`: see `forecast_all`
- `output_vars::Vector{Symbol}`: see `forecast_one`

#### Method 2:

- `meansbands_input_files::Dict{Symbol,S}`: dictionary mapping an output_var to the filename
  containing forecasts for that output_var (where `S<:AbstractString)`. Keys should be one of the following:
  `:histpseudo, :forecastpseudo, :shockdecpseudo, :forecastobs, :shockdecobs`.

### Keyword Arguments

- `density_bands::Vector{T}`: a vector of percent values (between 0 and 1) for which to compute density bands.

- `minimize::Bool`: if `true`, choose shortest interval, otherwise just chop off lowest and
  highest (percent/2)

- `forecast_string::S`: forecast identifier string (the value \"fcid=value\" in
  the forecast output filename). Required when `input_type == :subset`.

- `population_forecast_file::S`: if you have population forecast data,
  this is the filepath identifying where it is stored. In the method
  that accepts a model object, if `use_population_forecast(m) ==
  true`, the following file is used, if it exists:
  `inpath(m, \"data\", \"population_forecast_(data_vintage(m)).csv\")`

- `verbose`: level of error messages to be printed to screen. Options
  are `:none`, `:low`, `:high`

### Method 1:

- `load_dataset::Bool`: indicates whether or not to load the
  data using `load_data(m)`. Loading the dataset is required only when using
  transformations that convert values from log-levels to growth
  rates. Defaults to `true.`

- `load_population_data::Bool`: indicates whether or not to load the
  population growth rate data. This is required only when a series
  requires either the `:loglevelto4qpct_annualized_percapita` or
  `:loglevelto4qpct_annualized` transformation.

#### Method 2:

- `y0_indexes::Dict{Symbol,Int}`: A `Dict` storing the mapping of products to the index
  of the period prior to the first period for which that
  product is computed. This is used to compute growth rates of
  forecasted or counterfactual variables, such as the deterministic
  trend. `y0_indexes[:forecast]` should correspond to the last
  historical period; `y0_indexes[:dettrend]` should correspond
  to the last presample period. It is required for only those
  observables and pseudoobservables that employ the
  `loglevelto4qpct_annualized_percapita` transformation.

- `output_dir::S`: Directory in which to write means and bands files. Defaults to `""`.

- `population_data::DataFrame`: `DataFrame` containing the columns
   `:dlfiltered_population_recorded` and `date`.
   (`dlfiltered_population_recorded` refers to the log difference of
   the filtered population series, or the computed filtered population
   growth rate). Defaults to an empty `DataFrame` and is required in
   the same cases as the `load_population_data` argument in Method 1.

- `population_mnemonic::Nullable{Symbol}`: The name of the series holding the desired
  population series in `population_data`. Defaults to a `Nullable{Symbol}()`.

- `population_forecast_file::S`: Name of file in which to find population_data.

- `data::Matrix{T}`: pre-loaded `nobs x nperiods` matrix containing the transformed data matrix.
"""
function means_bands_all{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                           cond_type::Symbol, output_vars::Vector{Symbol};
                                           density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                           minimize::Bool = false,
                                           forecast_string = "",
                                           verbose::Symbol = :low)

    ## Step 0: Determine full set of output_vars necessary for plotting desired results
    #          Specifically, if output_vars contains shockdecs but not trend or deterministic trends,
    #          add those

    output_vars = add_requisite_output_vars(output_vars)


    ## Step 1: Population data

    # Parse population mnemonic
    tmp = parse_population_mnemonic(m)[1]
    population_mnemonic = tmp == Symbol() ? Nullable{Symbol}() : Nullable(tmp)

    # Get filenames for population data and forecast
    if !isnull(population_mnemonic)
        vint = data_vintage(m)
        population_data_file     = inpath(m, "data", "population_data_levels_$vint.csv")
        population_forecast_file = inpath(m, "data", "population_forecast_$vint.csv")
    else
        population_data_file = ""
        population_forecast_file = ""
    end


    ## Step 2: Load main dataset (required for some transformations),
    ##         specify which period to use as first t-1 period for computing
    ##         growth rates, and which period is the last historical period

    # Load dataset
    data = df_to_matrix(m, load_data(m, verbose = :none))

    # Specify the t-1 period for each product
    y0_indexes = Dict{Symbol,Int}()
    products   = unique(map(get_product, output_vars))
    for prod in products
        y0_indexes[prod] = get_y0_index(m, prod)
    end


    ## Step 3: Fetch the model string and the input and output directories we'll read from/write to

    tmp          = rawpath(m, "forecast", "foo")
    input_dir    = dirname(tmp)
    model_string = replace(basename(tmp), "foo_", "")
    output_dir   = dirname(workpath(m, "forecast", "foo"))


    ## Step 4: We have everything we need; appeal to model-object-agnostic function

    means_bands_all(input_type, cond_type, output_vars,
                    input_dir, output_dir;
                    density_bands = density_bands, minimize = minimize,
                    model_string  = model_string,
                    forecast_string = forecast_string,
                    population_mnemonic = population_mnemonic,
                    population_data_file = population_data_file,
                    population_forecast_file = population_forecast_file,
                    y0_indexes = y0_indexes, data = data,
                    verbose = verbose)
end

function means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                               cond_type::Symbol,
                                               output_vars::Vector{Symbol},
                                               input_dir::S,
                                               output_dir::S;
                                               density_bands::Vector{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               minimize::Bool = false,
                                               model_string = "",
                                               forecast_string = "",
                                               population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                                               population_data_file = "",
                                               population_forecast_file = "",
                                               y0_indexes::Dict{Symbol,Int} = Dict{Symbol,Int}(),
                                               data = Matrix{T}(),
                                               verbose::Symbol = :low)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println()
        info("Computing means and bands for input_type = $input_type, cond_type = $cond_type...")
        println("Start time: $(now())")
        println("Means and bands will be saved in $output_dir")
    end

    ## Preliminary step: if string arguments aren't of type S, force them to be
    forecast_string = S(forecast_string)
    model_string    = S(model_string)

    ## Step 1: Filter population history and forecast and compute growth rates

    dlfiltered_population_data, dlfiltered_population_forecast =
        load_population_growth(population_data_file, population_forecast_file,
                               population_mnemonic; verbose = verbose)

    ## Step 2: Set up filenames for MeansBands input and output files.
    # MeansBands output filenames are the same as forecast output filenames, but with an "mb" prefix.

    mb_input_files = get_meansbands_input_files(input_dir, input_type, cond_type,
                                                output_vars; model_string = model_string,
                                                forecast_string = forecast_string)
    mb_output_files = get_meansbands_output_files(output_dir, input_type, cond_type, output_vars,
                                                  model_string = model_string,
                                                  forecast_string = forecast_string)

    ## Step 3: Compute means and bands for each output variable, and write to a file.

    for output_var in output_vars

        # Which product are we dealing with? Need this to index out of y0_indexes.
        prod     = get_product(output_var)
        y0_index = isempty(y0_indexes) ? -1 : y0_indexes[prod]

        # compute means and bands object
        mb = means_bands(input_type, cond_type, output_var,
                         mb_input_files, density_bands = density_bands, minimize = minimize,
                         forecast_string = forecast_string,
                         population_data = dlfiltered_population_data,
                         population_forecast = dlfiltered_population_forecast,
                         population_mnemonic = Nullable(:population_growth),
                         y0_index = y0_index, data = data)

        # write to file
        filepath = mb_output_files[output_var]
        dirpath  = dirname(filepath)
        !isdir(dirpath) ? mkdir(dirpath) : nothing
        jldopen(filepath, "w") do file
            write(file, "mb", mb)
        end

        if VERBOSITY[verbose] >= VERBOSITY[:high]
            println(" * Wrote $(basename(filepath))")
        end
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("\nComputation of means and bands complete: $(now())")
    end
end

"""
```
means_bands{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                                 cond_type::Symbol,
                                                 output_var::Symbol,
                                                 meansbands_input_files::Dict{Symbol,S};
                                                 density_bands::Vector{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                                 minimize::Bool = false,
                                                 forecast_string::S = "",
                                                 population_data = DataFrame(),
                                                 population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                                                 population_forecast = DataFrame(),
                                                 y0_index::Int = -1,
                                                 data = Matrix{T}())
```

Computes means and bands for a single `output_var`.

### Input Arguments

All inputs are exactly the same as the second
`means_bands_all` method, except that `output_var` is a single
`Symbol` rather than `Array{Symbol}`.
"""
function means_bands{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                                          cond_type::Symbol,
                                                          output_var::Symbol,
                                                          meansbands_input_files::Dict{Symbol,S};
                                                          density_bands::Vector{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                                          minimize::Bool = false,
                                                          forecast_string::S = "",
                                                          population_data = DataFrame(),
                                                          population_forecast = DataFrame(),
                                                          population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                                                          y0_index::Int = -1,
                                                          data = Matrix{T}())

    # Return only one set of bands if we read in only one draw
    if input_type in [:init, :mode, :mean]
        density_bands = [.5]
    end

    ## Step 1: Determine the class of variable we are working with (pseudos? observables? etc)
    ##         and the product we are computing (forecast? history? shockdec?)

    class = get_class(output_var)
    product = get_product(output_var)

    ## Step 2: Read in raw forecast output and metadata (transformations, mappings from symbols to indices, etc)
    # open correct input file
    forecast_output_file = meansbands_input_files[output_var]
    fcast_output, metadata, mb_metadata, transforms =
        get_mb_metadata(input_type, cond_type, output_var, forecast_output_file;
                        forecast_string = forecast_string)

    date_list         = collect(keys(mb_metadata[:date_inds]))
    variable_indices  = mb_metadata[:indices]
    population_series = get_mb_population_series(product, population_mnemonic, population_data,
                                                 population_forecast, date_list)

    means, bands = if product in [:hist, :forecast, :dettrend, :trend, :forecast4q,
                                  :bddforecast, :bddforecast4q]

        # indicate whether we are getting a 4q forecast or regular annualized forecast
        fourquarter = product in [:forecast4q, :bddforecast4q]

        # The `fcast_output` for trends only is of size `ndraws` x `nvars`. We
        # need to use `repeat` below because population adjustments will be
        # different in each period. Now we have something of size `ndraws` x
        # `nvars` x `nperiods`
        if product == :trend
            fcast_output = repeat(fcast_output, outer = [1, 1, length(date_list)])
        end

        compute_means_bands(fcast_output, transforms, variable_indices;
                            date_list = date_list, data = data,
                            population_series = population_series, y0_index =
                            y0_index, density_bands = density_bands,
                            minimize = minimize, fourquarter = fourquarter)

    elseif product in [:shockdec, :irf]

        # get shock indices
        mb_metadata[:shock_indices] = metadata[:shock_indices]

        # compute means and bands for shock decomposition
        compute_means_bands(fcast_output, transforms, variable_indices,
                            metadata[:shock_indices]; date_list = date_list,
                            data = data, population_series = population_series,
                            y0_index = y0_index, density_bands = density_bands,
                            minimize = minimize)
    end

    return MeansBands(mb_metadata, means, bands)
end

"""
```
compute_means_bands(fcast_output::Array{T, 3}, transforms, variable_inds;
    date_list = [], data = [], population_series = [], y0_index = -1,
    density_bands = [0.5, 0.6, 0.7, 0.8, 0.9], minimize = false)
```

Compute means and bands for 3-dimensional (draw x variable x time) products:
histories, forecasts, deterministic trends, and trends.
"""
function compute_means_bands{T<:AbstractFloat}(fcast_output::Array{T, 3},
                                               transforms::Dict{Symbol,Symbol},
                                               variable_inds::Dict{Symbol,Int};
                                               date_list::Vector{Date} = Vector{Date}(),
                                               data::Matrix{T} = Matrix{T}(),
                                               population_series = Vector{T}(),
                                               y0_index::Int = -1,
                                               density_bands::Array{Float64} = [0.5,0.6,0.7,0.8,0.9],
                                               minimize::Bool = false,
                                               fourquarter::Bool = false)

    # Set up means and bands structures
    means = DataFrame(date = date_list)
    bands = Dict{Symbol,DataFrame}()

    # For each series (i.e. for each pseudoobs, obs, or state):
    # 1. Apply the appropriate transform
    # 2. Compute means and density bands of transformed output
    # 3. Add to DataFrames
    for (var, ind) in variable_inds

        # apply transformation to all draws
        transform = parse_transform(transforms[var])
        fcast_series = squeeze(fcast_output[:, ind, :], 2)

        transformed_fcast_output = if fourquarter
            transform4q = get_transform4q(transform)

            # transform
            result = if transform4q in [logtopct_4q_percapita]
                # we use y0_index+1 when we want to sum the last 4 periods
                hist_data = squeeze(data[ind, y0_index+1:end],1)
                transform4q(fcast_series, hist_data, population_series)
            elseif transform4q in [loglevelto4qpct_4q_percapita]
                # we use y0_index for computing growth rates
                hist_data = squeeze(data[ind, y0_index:end],1)
                transform4q(fcast_series, hist_data, population_series)
            elseif transform4q in [quartertoannual]
                transform4q(fcast_series)
            elseif transform4q in [logtopct_4q]
                # we use y0_index+1 when we want to sum the last 4 periods
                hist_data = squeeze(data[ind, y0_index+1:end], 1)
                transform4q(fcast_series, hist_data)
            elseif transform4q in [identity]
                fcast_series
            else
                error("Please provide an invocation for $transform4q in $(@__FILE__())")
            end

            result
        else
            result = if transform in [logtopct_annualized_percapita]
                transform(fcast_series, population_series)
            elseif transform in [loglevelto4qpct_annualized_percapita]
                hist_data = data[ind, y0_index]
                transform(fcast_series, hist_data, population_series)
            else
                transform(fcast_series)
            end

            result
        end

        # compute the mean and bands across draws and add to dataframe
        means[var] = vec(mean(transformed_fcast_output,1))
        bands[var] = find_density_bands(transformed_fcast_output, density_bands, minimize = minimize)
        bands[var][:date] = date_list
    end

    return means, bands
end


"""
```
compute_means_bands(fcast_output::Array{T, 4}, transforms, variable_inds, shock_inds;
    date_list = [], data = [], population_series = [], y0_index = -1,
    density_bands = [0.5, 0.6, 0.7, 0.8, 0.9], minimize = false)
```

Compute means and bands for 4-dimensional (draw x variable x time x shock)
products: shock decompositions and impulse responses. Note the extra input
argument `shock_inds` which is not present in the 2-dimensional method of
`compute_means_bands`.
"""
function compute_means_bands{T<:AbstractFloat}(fcast_output::Array{T, 4},
                                               transforms::Dict{Symbol,Symbol},
                                               variable_inds::Dict{Symbol,Int},
                                               shock_inds::Dict{Symbol,Int};
                                               date_list::Vector{Date} = Vector{Date}(),
                                               data::Matrix{T} = Matrix{T}(),
                                               population_series = Vector{T}(),
                                               y0_index::Int = -1,
                                               density_bands::Array{Float64} = [0.5,0.6,0.7,0.8,0.9],
                                               minimize::Bool = false)

    # Set up means and bands structures
    if !isempty(date_list)
        sort!(date_list)
        means = DataFrame(date = date_list)
    else
        means = DataFrame()
    end
    bands = Dict{Symbol,DataFrame}()

    # For each element of shock x variable (variable = each pseudoobs, obs, or state):
    # 1. Apply the appropriate transform
    # 2. Compute means and density bands of transformed output
    # 3. Add to DataFrames
    for (shock, shock_ind) in shock_inds
        for (var, var_ind) in variable_inds
            transform = parse_transform(transforms[var])
            fcast_series = squeeze(fcast_output[:, var_ind, :, shock_ind], 2)

            transformed_fcast_output = if transform in [logtopct_annualized_percapita]
                transform(fcast_series, population_series)
            elseif transform in [loglevelto4qpct_annualized_percapita]
                hist_data = data[var_ind, y0_index]
                transform(fcast_series, hist_data, population_series)
            else
                transform(fcast_series)
            end

            means[symbol("$var$DSGE_SHOCKDEC_DELIM$shock")] = vec(mean(transformed_fcast_output,1))

            bands_one = find_density_bands(transformed_fcast_output, density_bands; minimize = minimize)
            if !isempty(date_list)
                bands_one[:date] = date_list
            end
            bands[symbol("$var$DSGE_SHOCKDEC_DELIM$shock")] = bands_one
        end
    end

    return means, bands
end
