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
- `subset_string`: subset identifier string (the value "subs=value" in the forecast output file identifier string). Only
  to be used when `input_type == :subset`.
- `population_forecast_file:` if you have population forecast data,
  this is the filepath identifying where it is stored. In the method
  that accepts a model object, if `use_population_forecast(m) ==
  true`, the following file is used, if it exists:
"""
function compute_means_bands_all{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9],
                                               subset_string = "", population_forecast_file = "",
                                               load_dataset::Bool = true)

    # Get population forecast file
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

    # Load main dataset - required for some transformations
    data = load_dataset ? df_to_matrix(m, load_data(m)) : Matrix{T}()
    hist_end_index = if cond_type in [:semi, :full]
        n_mainsample_periods(m) - 1
    else
        n_mainsample_periods(m)
    end

    # Get names of files that the forecast wrote
    forecast_output_files = DSGE.get_output_files(m, "forecast", input_type,
                                                  output_vars, cond_type, subset_string = subset_string)


    compute_means_bands_all(input_type, output_vars, cond_type, forecast_output_files,
                            density_bands = density_bands, subset_string = subset_string,
                            output_dir = workpath(m,"forecast",""),
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
                                               population_forecast_file = "",
                                               hist_end_index::Int = 0,
                                               data = Matrix{T}())

    # read in population forecast file (if it's provided)
    pop_forecast = if !isempty(population_forecast_file)
        println("Reading population forecasts from $population_forecast_file\n")
        convert(Array{Float64}, readtable(population_forecast_file)[:POPULATION])
    else
        warn("No population forecast supplied.")
        []
    end

    # set mb_output_vars so mb is the prefix to all meansbands output files
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

    for (i,output_var) in enumerate(output_vars)

        # compute means and bands object
        mb = compute_means_bands(input_type, output_var, cond_type,
                                 forecast_output_files, density_bands = density_bands,
                                 subset_string = subset_string,
                                 population_forecast = pop_forecast,
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
                                                    population_forecast = [],
                                                    hist_end_index::Int = 0,
                                                    data = Matrix())

    # Determine whether we should compute means and bands for pseudos or observables
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

    # TODO
    if product == :shockdec
        return compute_means_bands_shockdec(fcast_output, transforms, inds, date_list, date_inds_order)
    end

    # Ensure population forecast is same length as fcast_output.
    # For forecasts, the third dimension of the fcast_output matrix is the number of periods.
    population_forecast = if product in [:forecast]
        nperiods = size(fcast_output, 3)
        resize_population_forecast(population_forecast, nperiods)
    else
        population_forecast
    end

    # make DataFrames for means and bands
    sort!(date_list)
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
            println("size of data: $(size(data))")
            println(size(squeeze(fcast_output[:,ind,date_inds_order],2)))
            println("end_index: $hist_end_index")
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


function compute_means_bands_shockdec(fcast_output, transforms, inds, date_list, date_inds_order)

    # do *not* apply transform

    # deal with pop forecast

    # include shock_inds in mb.metadata
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