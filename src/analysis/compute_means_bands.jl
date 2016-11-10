"""
```
compute_means_bands_all{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
    output_vars::Vector{Symbol}, cond_type::Symbol;
    density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9], subset_string = "",
    population_forecast_file = "", verbose::Symbol = :low)

compute_means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
    output_vars::Vector{Symbol}, cond_type::Symbol, forecast_files::Dict{Symbol,S};
    density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9], subset_string = "",
    output_dir = "", population_forecast_file = "", verbose::Symbol = :low)

```

Computes means and bands for pseudoobservables and observables, and writes
results to a file. Two methods are provided. The method that accepts a model
object as an argument uses the model's settings to infer `forecast_files`; then
appeals to the second method. Users can optionally skip construction of a model
object and manually enter `forecast_files`.

### Input Arguments

- `m`: model object
- `input_type`: see `forecast_all`
- `output_vars`: see `forecast_one`
- `cond_type`: see `forecast_all`
- `forecast_files`: dictionary mapping an output_var to the filename
  containing forecasts for that output_var. Keys should be one of the following:
  `:histpseudo, :forecastpseudo, :shockdecpseudo, :forecastobs, :shockdecobs`.

### Keyword Arguments

- `density_bands`: a vector of percent values (between 0 and 1) for which to
  compute density bands.
- `subset_string`: subset identifier string (the value "subs=value" in the
  forecast output file identifier string). Only to be used when
  `input_type == :subset`.
- `population_forecast_file:` if you have population forecast data, this is the
  filepath identifying where it is stored. In the method that accepts a model
  object, if `use_population_forecast(m) == true`, the following file is used,
  if it exists:
- `verbose::Symbol`: desired frequency of function progress messages printed to
  standard out. One of `:none`, `:low`, or `:high`
"""
function compute_means_bands_all{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               subset_string = "", load_dataset::Bool = true,
                                               load_population_data::Bool = true,
                                               population_forecast_file = "",
                                               verbose::Symbol = :low)

    ## Step 1: Get population forecast file
    if isempty(population_forecast_file)
        population_forecast_file = if use_population_forecast(m)
            inpath(m, "data", "population_forecast_$(data_vintage(m)).csv")
        else
            ""
        end
    end

    # load population level data, which was saved in load_data_levels
    if load_population_data
        level_data = read_population_data(m)
    end

    population_mnemonic = Nullable(parse_population_mnemonic(m)[1])

    ## Step 2: Load main dataset - required for some transformations
    data = load_dataset ? df_to_matrix(m, load_data(m; verbose = :none)) : Matrix{T}()

    ## Step 3: Get names of files that the forecast wrote
    forecast_output_files = DSGE.get_output_files(m, "forecast", input_type,
                                                  output_vars, cond_type, subset_string = subset_string)

    ## Step 4: We have everything we need; appeal to model-object-agnostic function
    compute_means_bands_all(input_type, output_vars, cond_type, forecast_output_files,
                            density_bands = density_bands, subset_string = subset_string,
                            output_dir = workpath(m,"forecast",""),
                            population_data = level_data,
                            population_mnemonic = population_mnemonic,
                            population_forecast_file = population_forecast_file,
                            data = data, verbose = verbose)
end

function compute_means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                               output_vars::Vector{Symbol},
                                               cond_type::Symbol,
                                               forecast_output_files::Dict{Symbol,S};
                                               density_bands::Vector{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               subset_string = "",
                                               output_dir = "",
                                               population_data::DataFrame = DataFrame(),
                                               population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                                               population_forecast_file = "",
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
            isempty(population_data) ? warn("No population data provided") : nothing
            isnull(population_mnemonic) ? warn("No population mnemonic provided") : nothing

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

        # compute means and bands object
        mb = compute_means_bands(input_type, output_var, cond_type,
                                 forecast_output_files, density_bands = density_bands,
                                 subset_string = subset_string,
                                 population_data = dlfiltered_population_data,
                                 population_mnemonic = Nullable(:population_growth),
                                 population_forecast = dlfiltered_population_forecast,
                                 hist_end_index = hist_end_index,
                                 data = data)

        # write to file
        filepath = mb_files[output_var]
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

function compute_means_bands{S<:AbstractString}(input_type::Symbol,
                                                output_var::Symbol,
                                                cond_type::Symbol,
                                                forecast_output_files::Dict{Symbol,S};
                                                density_bands = [0.5, 0.6, 0.7, 0.8, 0.9],
                                                subset_string::S = "",
                                                population_data = DataFrame(),
                                                population_mnemonic::Nullable{Symbol} = Nullable{Symbol}(),
                                                population_forecast = DataFrame(),
                                                hist_end_index::Int = 0,
                                                data = Matrix())

    # Return only one set of bands if we read in only one draw
    if input_type in [:init, :mode, :mean]
        density_bands = [.5]
    end

    ## Step 1: Determine the class of variable we are working with (pseudos? observables? etc)
    ##         and the product we are computing (forecast? history? shockdec?)
    class = if contains(string(output_var), "pseudo")
        :pseudo
    elseif contains(string(output_var), "obs")
        :obs
    elseif contains(string(output_var), "state")
        :state
    elseif contains(string(output_var), "shock")
        :shock
    end

    product = if contains(string(output_var), "hist")
        :hist
    elseif contains(string(output_var), "forecast")
        :forecast
    elseif contains(string(output_var), "shockdec")
        :shockdec
    end

    ## Step 2: Read in raw forecast output and metadata (transformations, mappings from symbols to indices, etc)
    # open correct input file
    forecast_output_file = forecast_output_files[output_var]
    metadata, fcast_output = jldopen(forecast_output_file, "r") do jld
        read_forecast_metadata(jld), DSGE.read_darray(jld)
    end

    transforms, variable_indices, date_indices = if class == :pseudo
        metadata[:pseudoobservable_revtransforms], metadata[:pseudoobservable_indices],
        metadata[:date_indices]
    elseif class == :obs
        metadata[:observable_revtransforms], metadata[:observable_indices], metadata[:date_indices]
    else
        error("means and bands are only calculated for observables and pseudo-observables")
    end

    # make sure date lists are valid
    date_list          = collect(keys(date_indices))   # unsorted array of actual dates
    date_indices_order = collect(values(date_indices)) # unsorted array of date indices
    check_consistent_order(date_list, date_indices_order)
    sort!(date_list, by = x -> date_indices[x])
    sort!(date_indices_order)

    # get population mnemonic
    mnemonic = isnull(population_mnemonic) ? Symbol() : get(population_mnemonic)

    # Ensure population forecast is same length as fcast_output.
    # For forecasts, the third dimension of the fcast_output matrix is the number of periods.
    population_forecast = if product in [:forecast, :shockdec]
        n_fcast_periods = size(fcast_output, 3)
        resize_population_forecast(population_forecast, n_fcast_periods,
                                   population_mnemonic = mnemonic)
    end

    mb_metadata = Dict{Symbol,Any}(
                   :para       => input_type,
                   :cond_type  => cond_type,
                   :product    => product,
                   :class      => class,
                   :indices    => variable_indices,
                   :subset_string => subset_string)

    means, bands = if product in [:shockdec]
        # make sure population series corresponds with saved shockdec dates
        shockdec_start = date_list[1]
        shockdec_end   = date_list[end]

        start_ind = find(population_data[:date] .== shockdec_start)[1]
        end_ind   = find(population_forecast[:date] .== shockdec_end)[1]

        # concatenate population histories and forecasts together
        population_series = if isempty(end_ind)
            convert(Vector{Float64}, population_data[start_ind:end, mnemonic])
        else
            tmp = [population_data[start_ind:end, mnemonic]; population_forecast[1:end_ind, mnemonic]]
            convert(Vector{Float64}, tmp)
        end

        # get shock indices
        mb_metadata[:shock_indices] = metadata[:shock_indices]

        # compute means and bands for shock decomposition
        compute_means_bands_shockdec(fcast_output, transforms,
                                     variable_indices, metadata[:shock_indices], date_list,
                                     data = data, population_series = population_series,
                                     density_bands = density_bands)
    else
        # make DataFrames for means and bands
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        # for each series (ie each pseudoobs, each obs, or each state):
        # 1. apply the appropriate transform
        # 2. add to DataFrame
        for (series, ind) in variable_indices
            # apply transformation to all draws
            transform = parse_transform(transforms[series])
            fcast_series = squeeze(fcast_output[:, ind, :], 2)
            transformed_fcast_output = if transform in [logtopct_annualized_percapita]
                pop_fcast = convert(Vector{Float64}, population_forecast[mnemonic])
                transform(fcast_series, pop_fcast)
            elseif transform in [loglevelto4qpct_annualized_percapita]
                pop_fcast = convert(Vector{Float64}, population_forecast[mnemonic])
                hist_series = vec(data[ind, :])
                transform(fcast_series, hist_series, pop_fcast)
            else
                map(transform, fcast_series)
            end

            # compute the mean and bands across draws and add to dataframe
            means[series] = vec(mean(transformed_fcast_output,1))
            bands[series] = find_density_bands(transformed_fcast_output, density_bands, minimize=false)
            bands[series][:date] = date_list
        end

        means, bands
    end

    return MeansBands(mb_metadata, means, bands)
end

