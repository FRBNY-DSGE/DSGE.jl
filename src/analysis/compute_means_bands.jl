"""
```
compute_means_bands_all{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               subset_string = "", population_forecast_file = "")

compute_means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol, output_vars::Vector{Symbol},
                                               cond_type::Symbol, forecast_files::Dict{Symbol,S};
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               subset_string = "", output_dir = "",
                                               population_forecast_file = "")

```

Computes means and bands for pseudoobservables and observables, and writes results to a file. Two
methods are provided. The method that accepts a model object as an
argument uses the model's settings to infer `forecast_files`; then
appeals to the second method. Users can optionally skip construction
of a model object and manually enter `forecast_files`.

### Input Arguments

- `m`: model object
- `input_type`: see `forecast_all`
- `output_vars`: see `forecast_one`
- `cond_type`: see `forecast_all`
- `forecast_files`: dictionary mapping an output_var to the filename
  containing forecasts for that output_var. Keys should be one of the following:
  `:histpseudo, :forecastpseudo, :shockdecpseudo, :forecastobs, :shockdecobs`.

### Keyword Arguments

- `density_bands`: a vector of percent values (between 0 and 1) for which to compute density bands.
- `subset_string`: subset identifier string (the value "subs=value" in
  the forecast output file identifier string). Only to be used when
  `input_type == :subset`.
- `population_forecast_file:` if you have population forecast data,
  this is the filepath identifying where it is stored. In the method
  that accepts a model object, if `use_population_forecast(m) ==
  true`, the following file is used, if it exists:
"""
function compute_means_bands_all{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               subset_string = "", load_dataset::Bool = true
                                               load_population_data::Bool = false,
                                               population_forecast_file = "",
                                               verbose::Symbol = :low)

    ## Step 1: Get population forecast file
    population_forecast_file = if isempty(population_forecast_file)
        guess = if use_population_forecast(m)
            inpath(m, "data", "population_forecast_$(data_vintage(m)).csv")
        else
            ""
        end

        guess
    else
        population_forecast_file
    end

    # load population level data, which happens in load_data_levels
    level_data = if load_population_data
        load_data_levels(m, verbose = verbose)
    end

    ## Step 2: Load main dataset - required for some transformations
    data = load_dataset ? df_to_matrix(m, load_data(m)) : Matrix{T}()
    hist_end_index = if cond_type in [:semi, :full]
        n_mainsample_periods(m) - 1
    else
        n_mainsample_periods(m)
    end

    ## Step 3: Get names of files that the forecast wrote
    forecast_output_files = DSGE.get_output_files(m, "forecast", input_type,
                                                  output_vars, cond_type, subset_string = subset_string)

    ## Step 4: We have everything we need; appeal to model-object-agnostic function
    compute_means_bands_all(input_type, output_vars, cond_type, forecast_output_files,
                            density_bands = density_bands, subset_string = subset_string,
                            output_dir = workpath(m,"forecast",""),
                            population_data = level_data,
                            population_forecast_file = population_forecast_file,
                            hist_end_index = hist_end_index, data = data)
end
function compute_means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                               output_vars::Vector{Symbol},
                                               cond_type::Symbol,
                                               forecast_output_files::Dict{Symbol,S};
                                               density_bands::Vector{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               subset_string = "",
                                               output_dir = "",
                                               population_data::DataFrame = DataFrame(),
                                               population_mnemonic::Symbol = symbol(""),
                                               population_forecast_file = "",
                                               hist_end_index::Int = 0,
                                               data = Matrix{T}())


    ## Step 1: Filter population history and forecast and compute growth rates

    dlfiltered_population_data, dlfiltered_population_forecast =
        if !isempty(population_data) && !isempty(population_mnemonic)

            # get all of the population data
            population_data = transform_population_data(population_data, population_mnemonic,
                                                        population_forecast_file = population_forecast_file)

            population_data[:,[:date, :dlfiltered_population_recorded]],
            population_data[:,[:date, :dlfiltered_population_forecast]]
        else
            DataFrame(), DataFrame()
        end

    ## Step 2: Set up filenames for MeansBands output files.
    # MeansBands output filenames are the same as forecast output filenames, but with an "mb" prefix.
    mb_output_vars = [symbol("mb$x") for x in output_vars]

    mb_files = Dict{Symbol,AbstractString}()
    for x in keys(forecast_output_files)
        fn = forecast_output_files[x]
        base = "mb" * basename(fn)
        mb_files[x] = if isempty(output_dir)
            dir  = dirname(fn)
            joinpath(dir,base)
        else
            joinpath(output_dir,base)
        end
    end

    ## Step 3: Compute means and bands for each output variable, and write to a file.
    for (i,output_var) in enumerate(output_vars)

        # compute means and bands object
        mb = compute_means_bands(input_type, output_var, cond_type,
                                 forecast_output_files, density_bands = density_bands,
                                 subset_string = subset_string,
                                 population_data = dlfiltered_population_data,
                                 population_forecast = dlfiltered_population_forecast
                                 hist_end_index = hist_end_index,
                                 data = data)

        # write to file
        filepath = mb_files[output_vars[i]]
        jldopen(filepath, "w") do file
            write(file, "mb", mb)
        end
        println("wrote means and bands for ($output_var, $cond_type) to $filepath.\n")
    end
end

function compute_means_bands{S<:AbstractString}(input_type::Symbol,
                                                    output_var::Symbol,
                                                    cond_type::Symbol,
                                                    forecast_output_files::Dict{Symbol,S};
                                                    density_bands = [0.5, 0.6, 0.7, 0.8, 0.9],
                                                    subset_string::S = "",
                                                    population_data = DataFrame(),
                                                    population_forecast = DataFrame(),
                                                    hist_end_index::Int = 0,
                                                    data = Matrix())

    # Return only one set of bands if we read in only one draw
    if input_type in [:init, :mode, :mean]
        density_bands = [1.]
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

    print("* Computing means and bands for $output_var...")

    ## Step 2: Read in raw forecast output and metadata (transformations, mappings from symbols to indices, etc)
    # open correct input file
    forecast_output_file = forecast_output_files[output_var]
    metadata, fcast_output = jldopen(forecast_output_file, "r") do jld

        # read metadata
        metadata = read_forecast_metadata(jld)

        # read the DArray using read_darray
        fcast_output = DSGE.read_darray(jld)

        metadata, fcast_output
    end

    transforms, inds, date_inds = if class == :pseudo
        metadata[:pseudoobservable_revtransforms], metadata[:pseudoobservable_indices],
        metadata[:date_indices]
    elseif class == :obs
        metadata[:observable_revtransforms], metadata[:observable_indices], metadata[:date_indices]
    else
        error("means and bands are only calculated for observables and pseudoobservables")
    end

    # make sure date lists are valid
    date_list       = collect(keys(date_inds))   # unsorted array of actual dates
    date_inds_order = collect(values(date_inds)) # unsorted array of date indices
    check_consistent_order(date_list, date_inds_order)
    sort!(date_list, by = x -> date_inds[x])
    sort!(date_inds_order)

    # Ensure population forecast is same length as fcast_output.
    # For forecasts, the third dimension of the fcast_output matrix is the number of periods.
    population_forecast = if product in [:forecast]
        nperiods = size(fcast_output, 3)
        resize_population_forecast(population_forecast, nperiods)
    else
        population_forecast
    end

    if product == :shockdec

        # make sure population series corresponds with saved shockdec dates
        shockdec_start = date_list[1]
        shockdec_end   = date_list[end]

        start_ind = find(population_data[:date] .== shockdec_start)
        end_ind   = find(population_forecast[:date] .== shockdec_end)

        # concatenate population histories and forecasts together
        population_series = if isempty(end_ind)
            population_data[start_ind:end]
        else
            [population_data[start_ind:end]; population_forecast[1:end_ind]]
        end

        # compute means and bands for shock decomposition
        return compute_means_bands_shockdec(fcast_output[:,:,date_inds_order,:], transforms, inds, date_list,
                                            data = data, population_series = population_series,
                                            hist_end_index = hist_end_index)
    end

    # make DataFrames for means and bands
    means = DataFrame(date = date_list)
    bands = Dict{Symbol,DataFrame}()

    # for each series (ie each pseudoobs, each obs, or each state):
    # 1. apply the appropriate transform
    # 2. add to DataFrame
    for (series, ind) in inds
        # apply transformation to all draws
        transform = parse_transform(transforms[series])
        ex = if transform in [:logtopct_annualized]
            Expr(:call, transform, squeeze(fcast_output[:,ind,date_inds_order],2), population_forecast')
        elseif transform in [:loglevelto4qpct_annualized]
            Expr(:call, transform, squeeze(fcast_output[:,ind,date_inds_order],2), data[ind,:], hist_end_index)
        else
            Expr(:call, :map, transform, squeeze(fcast_output[:,ind,date_inds_order],2))
        end

        transformed_fcast_output = eval(ex)

        # compute bands
        bands_one = find_density_bands(transformed_fcast_output, density_bands, minimize=false)
        bands_one[:date] = date_list

        # compute the mean and bands across draws and add to dataframe
        means[series] = vec(mean(transformed_fcast_output,1))
        bands[series] = bands_one
    end

    mb_metadata = Dict{Symbol,Any}(
                   :para       => input_type,
                   :cond_type  => cond_type,
                   :product    => product,
                   :class      => class,
                   :indices    => inds,
                   :subset_string => subset_string)

    return MeansBands(mb_metadata, means, bands)
end

"""
```
compute_means_bands_shockdec(fcast_output, transforms, var_inds, shock_inds, date_list,
                             data = [], population_forecast = [], hist_end_index = 0)
```


### Inputs

- `fcast_output`: an `ndraws` x `nvars` x `nperiods` x `nshocks`
  matrix of shock decompositions from the output of the forecast

"""

function compute_means_bands_shockdec(fcast_output, transforms, var_inds, shock_inds, date_list;
                                      data = [], population_data = [], population_forecast = [], hist_end_index = 0)

    # set up means and bands structures
    sort!(date_list)
    means = DataFrame(date = date_list)
    bands = Dict{Symbol,DataFrame}()

    # for each element of shock x variable (variable = each pseudoobs, obs, or state):
    # 1. apply the appropriate transform
    # 2. add to DataFrame
    for (shock, shock_ind) in shock_inds

        for (var, var_ind) in var_inds

            for period in [:past, :future]
                transform = parse_transform(transforms[var])
                ex = if transform in [:logtopct_annualized]
                    Expr(:call, transform, squeeze(fcast_output[:,var_ind,:,shock_ind],2), population_series')
                elseif transform in [:loglevelto4qpct_annualized]
                    Expr(:call, transform, squeeze(fcast_output[:,var_ind,:,shock_ind],2),
                         data[var_ind,:], hist_end_index)
                else
                    Expr(:call, :map, transform, squeeze(fcast_output[:,ind,date_inds_order],2))
                end
            end

        end

    end
    # apply transform to historical and forecast periods separately usind different dlpop values

    # include shock_inds in mb.metadata
    mb.metadata[:shock_inds] = shock_inds

    error("Todo: implement compute_means_bands_shockdec")
end

"""
```
read_forecast_metadata(file::JLD.JldFile)
```

Read metadata from forecast output files. This includes dictionaries mapping dates, as well as state, observable,
pseudo-observable, and shock names, to their respective indices in the saved
forecast output array. The saved dictionaries include:

- `date_indices::Dict{Date, Int}`: saved for all forecast outputs
- `state_names::Dict{Symbol, Int}`: saved for `var in [:histstates, :forecaststates, :shockdecstates]`
- `observable_names::Dict{Symbol, Int}`: saved for `var in [:forecastobs, :shockdecobs]`
- `observable_revtransforms::Dict{Symbol, Symbol}`: saved identifiers for reverse transforms used for observables
- `pseudoobservable_names::Dict{Symbol, Int}`: saved for `var in [:histpseudo, :forecastpseudo, :shockdecpseudo]`
- `pseudoobservable_revtransforms::Dict{Symbol, Symbol}`: saved identifiers for reverse transforms used for pseudoobservables
- `shock_names::Dict{Symbol, Int}`: saved for `var in [:histshocks, :forecastshocks, :shockdecstates, :shockdecobs, :shockdecpseudo]`
"""
function read_forecast_metadata(file::JLD.JldFile)
    metadata = Dict{Symbol, Any}()
    for field in names(file)
        metadata[symbol(field)] = read(file, field)
    end

    metadata
end

"""
```
parse_transform(t::Symbol)
```

Parse the module name out of a Symbol to recover the transform associated with an observable or pseudoobservable.
"""
function parse_transform(t::Symbol)
    symbol(split(string(t),".")[end])
end

"""
```
check_consistent_order(l1, l2)
```

Checks to make sure that l1 and l2 are ordered consistently.
"""
function check_consistent_order(l1, l2)
    @assert length(l1) == length(l2)

    # cache old pairs in a Dict
    original_pairs = Dict{eltype(l1), eltype(l2)}()
    for (i,item) in enumerate(l1)
        original_pairs[item] = l2[i]
    end

    # sort
    l1_sorted = sort(l1)
    l2_sorted = sort(l2)

    # make sure each pair has same pair as before
    for i in length(l1)
        @assert original_pairs[l1_sorted[i]] == l2_sorted[i] "Lists not consistently ordered"
    end

    return true
end