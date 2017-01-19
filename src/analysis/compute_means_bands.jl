"""
```
means_bands_all(m, input_type, cond_type, output_vars;
                density_bands = [0.5, 0.6, 0.7, 0.8, 0.9], minimize = false,
                forecast_string = "", load_dataset = true, load_population_data = true,
                population_forecast_file = "", verbose :low)


means_bands_all(input_type, cond_type, output_vars, forecast_output_files;
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

- `forecast_output_files::Dict{Symbol,S}`: dictionary mapping an output_var to the filename
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
                                           forecast_string = "", load_dataset::Bool = true,
                                           load_population_data::Bool = true,
                                           population_forecast_file = "",
                                           verbose::Symbol = :low)

    ## Step 0: Determine full set of output_vars necessary for plotting desired results
    #          Specifically, if output_vars contains shockdecs but not trend or deterministic trends,
    #          add those

    output_vars = add_requisite_output_vars(output_vars)

    ## Step 1: Load population data in levels

    # get population forecast file
    if isempty(population_forecast_file)
        population_forecast_file = if use_population_forecast(m)
            inpath(m, "data", "population_forecast_$(data_vintage(m)).csv")
        else
            ""
        end
    end

    # load population level data, which was saved in load_data_levels
    level_data = if load_population_data
        read_population_data(m)
    else
        DataFrame()
    end

    # reformat population_mnemonic
    population_mnemonic = Nullable(parse_population_mnemonic(m)[1])

    ## Step 2: Load main dataset (required for some transformations),
    ##         specify which period to use as first t-1 period for computing
    ##         growth rates, and which period is the last historical period

    data, y0_indexes = if load_dataset

        # load dataset
        data = df_to_matrix(m, load_data(m, verbose = :none))

        # specify the t-1 period for each product
        products = unique(map(get_product, output_vars))

        y0_indexes = Dict{Symbol,Int}()
        for prod in intersect(products, [:forecast])
            y0_indexes[prod] = index_forecast_start(m) - 1
        end
        for prod in intersect(products, [:shockdec, :dettrend, :trend])
            y0_indexes[prod] = index_shockdec_start(m) - 1
        end
        for prod in intersect(products, [:hist])
            y0_indexes[prod] = index_mainsample_start(m) - 1
        end

        data, y0_indexes
    else
        Matrix{T}(), Dict{Symbol,Nullable{Int}}()
    end

    ## Step 3: Get names of files that the forecast wrote
    forecast_output_files = DSGE.get_output_files(m, "forecast", input_type, cond_type,
                                                  output_vars; forecast_string = forecast_string)

    ## Step 4: We have everything we need; appeal to model-object-agnostic function
    means_bands_all(input_type, cond_type, output_vars, forecast_output_files,
                    density_bands = density_bands, minimize = minimize,
                    forecast_string = forecast_string,
                    output_dir = workpath(m,"forecast",""),
                    population_data = level_data,
                    population_mnemonic = population_mnemonic,
                    population_forecast_file = population_forecast_file,
                    y0_indexes = y0_indexes, data = data,
                    verbose = verbose)
end

function means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                               cond_type::Symbol,
                                               output_vars::Vector{Symbol},
                                               forecast_output_files::Dict{Symbol,S};
                                               density_bands::Vector{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               minimize::Bool = false,
                                               forecast_string = "",
                                               output_dir = "",
                                               population_data::DataFrame = DataFrame(),
                                               population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
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

    ## Step 1: Filter population history and forecast and compute growth rates

    dlfiltered_population_data, dlfiltered_population_forecast =
        if !(isempty(population_data) || isnull(population_mnemonic))
            # get all of the population data
            population_data, population_forecast =
                transform_population_data(population_data, get(population_mnemonic),
                                          population_forecast_file = population_forecast_file,
                                          verbose = :none)

            DataFrame(date = @data(convert(Array{Date}, population_data[:date])),
                                 population_growth = @data(convert(Array{Float64},
                                 population_data[:dlfiltered_population_recorded]))),

            DataFrame(date = @data(convert(Array{Date}, population_forecast[:date])),
                              population_growth = @data(convert(Array{Float64},
                              population_forecast[:dlfiltered_population_forecast])))

        else
            isempty(population_data) && VERBOSITY[verbose] >= VERBOSITY[:low] ?
                warn("No population data provided") : nothing

            isnull(population_mnemonic) && VERBOSITY[verbose] >= VERBOSITY[:low] ?
                warn("No population mnemonic provided") : nothing

            DataFrame(), DataFrame()
        end

    ## Step 2: Set up filenames for MeansBands output files.
    # MeansBands output filenames are the same as forecast output filenames, but with an "mb" prefix.
    mb_output_vars = [symbol("mb$x") for x in output_vars]

    mb_files = Dict{Symbol,AbstractString}()
    for (x, fn) in forecast_output_files
        base = "mb" * basename(fn)
        mb_files[x] = if isempty(output_dir)
            dir  = dirname(fn)
            joinpath(dir,base)
        else
            joinpath(output_dir,base)
        end
    end

    ## Step 3: Compute means and bands for each output variable, and write to a file.

    for output_var in output_vars

        # Which product are we dealing with? Need this to index out of y0_indexes.
        prod = get_product(output_var)
        if !isempty(y0_indexes)
            y0_index = y0_indexes[prod]
        else
            y0_index = -1
        end

        # compute means and bands object
        mb = means_bands(input_type, cond_type, output_var,
                         forecast_output_files, density_bands = density_bands,
                         minimize = minimize,
                         forecast_string = forecast_string,
                         population_data = dlfiltered_population_data,
                         population_mnemonic = Nullable(:population_growth),
                         population_forecast = dlfiltered_population_forecast,
                         y0_index = y0_index,
                         data = data)

        # write to file
        filepath = mb_files[output_var]
        dirpath  = dirname(filepath)
        if !isdir(dirpath)
            mkdir(dirpath)
        end
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
                                                 forecast_output_files::Dict{Symbol,S};
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
                                                          forecast_output_files::Dict{Symbol,S};
                                                          density_bands::Vector{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                                          minimize::Bool = false,
                                                          forecast_string::S = "",
                                                          population_data = DataFrame(),
                                                          population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                                                          population_forecast = DataFrame(),
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
    forecast_output_file = forecast_output_files[output_var]
    metadata, fcast_output = jldopen(forecast_output_file, "r") do jld
        read_forecast_metadata(jld), DSGE.read_darray(jld)
    end

    if class == :pseudo
        transforms       = metadata[:pseudoobservable_revtransforms]
        variable_indices = metadata[:pseudoobservable_indices]
    elseif class == :obs
        transforms       = metadata[:observable_revtransforms]
        variable_indices = metadata[:observable_indices]
    else
        error("Means and bands are only calculated for observables and pseudo-observables")
    end
    date_indices         = product == :irf ? Dict{Date,Int}() : metadata[:date_indices]

    # Make sure date lists are valid. This is vacuously true for trend and IRFs,
    # which are not time-dependent and hence have empty `date_indices`.
    date_list          = collect(keys(date_indices))   # unsorted array of actual dates
    date_indices_order = collect(values(date_indices)) # unsorted array of date indices
    check_consistent_order(date_list, date_indices_order)
    sort!(date_list, by = x -> date_indices[x])
    sort!(date_indices_order)

    # get population mnemonic
    mnemonic = isnull(population_mnemonic) ? Symbol() : get(population_mnemonic)

    # Ensure population series is same length as fcast_output.
    population_series = if product in [:forecast]

        # For forecasts, the third dimension of the fcast_output
        # matrix is the number of periods.

        n_fcast_periods = size(fcast_output, 3)
        population_series = resize_population_forecast(population_forecast, n_fcast_periods,
                                                       population_mnemonic = mnemonic)

        convert(Vector{Float64}, population_series[mnemonic])

    elseif product in [:shockdec, :dettrend, :trend]

        # For shockdecs, deterministic trend, and trend, we want to
        # make sure population series corresponds with the saved dates.

        start_date = date_list[1]
        end_date   = date_list[end]
        start_ind  = find(population_data[:date] .== start_date)[1]
        population_data = population_data[start_ind:end, mnemonic]

        # calculate number of periods that are in the future
        n_fcast_periods = length(date_list) - length(population_data)

        # Extend population forecast by the right number of periods
        population_forecast = resize_population_forecast(population_forecast, n_fcast_periods,
                                                         population_mnemonic = mnemonic)

        end_ind   = find(population_forecast[:date] .== end_date)[1]

        # concatenate population histories and forecasts together
        population_series = if isempty(end_ind)
            convert(Vector{Float64}, population_data)
        else
            tmp = [population_data; population_forecast[1:end_ind, mnemonic]]
            convert(Vector{Float64}, tmp)
        end

        population_series

    elseif product in [:hist]

        # For history, the population series is just the data

        convert(Vector{Float64}, population_data[mnemonic])
    end

    mb_metadata = Dict{Symbol,Any}(
                   :para       => input_type,
                   :cond_type  => cond_type,
                   :product    => product,
                   :class      => class,
                   :indices    => variable_indices,
                   :forecast_string => forecast_string,
                   :date_inds  => date_indices)

    means, bands = if product in [:hist, :forecast, :dettrend, :trend]

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
                            minimize = minimize)

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
                                               minimize::Bool = false)

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

        transformed_fcast_output = if transform in [logtopct_annualized_percapita]
            transform(fcast_series, population_series)
        elseif transform in [loglevelto4qpct_annualized_percapita]
            hist_data = data[ind, y0_index]
            transform(fcast_series, hist_data, population_series)
        else
            transform(fcast_series)
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
