"""
```
type MeansBands
```

Stores the means and bands associated with a class of vairables (observables, pseudoobservables, and states).

### Fields
- `key::Symbol`: indicates class of variable that this `MeansBands`
  object summarizes. Should be one of the following: `:pseudo`,
  `:obs`, `:state`
- `means::DataFrame`: DataFrame storing the mean time series for each
  variable in the class (e.g. one column for each element of
  `m.observables` if `key==:obs`, or one column for each
  pseudoobservable)
- `bands::Dict{Symbol,DataFrame}`: `bands` should have one key for
  each variable in the class. The corresponding DataFrame contains
  confidence bands for each variable. See `find_density_bands` for
  more information.
"""
type MeansBands
    key::Symbol
    means::DataFrame
    bands::Dict{Symbol,DataFrame}
end


"""
```
compute_means_bands{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])

compute_means_bands{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol, output_vars::Vector{Symbol},
                                               cond_type::Symbol, forecast_files::Dict{Symbol,S};
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])
```

Computes means and bands for pseudoobservables and observables. Two
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

### Outputs

`compute_means_bands` returns 2 `MeansBands` objects; one for observables and one for pseudoobservables.
"""
function compute_means_bands{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])

    # Get names of files that the forecast wrote
    forecast_output_files =  DSGE.get_output_files(m, input_type, output_vars, cond_type)
    compute_means_bands(input_type, output_vars, cond_type, forecast_output_files, density_bands = density_bands)
end
function compute_means_bands{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol, output_vars::Vector{Symbol},
                                               cond_type::Symbol, forecast_files::Dict{Symbol,S};
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])

    # make output dictionaries
    meansobs    = Dict{Symbol, DataFrame}()
    meanspseudo = Dict{Symbol, DataFrame}()
    bandsobs    = Dict{Symbol, Dict{Symbol,DataFrame}}()
    bandspseudo = Dict{Symbol, Dict{Symbol,DataFrame}}()

    # fill output dictionaries
    for output_var in output_vars  # histstates, forecastpseudo, shockdecobs, etc

        # Determine whether we should compute means and bands for pseudos or observables
        pseudo_flag = contains(string(output_var), "pseudo")
        obs_flag    = contains(string(output_var), "obs")

        # open correct input file
        forecast_output_file = forecast_output_files[output_var]
        metadata, data = jldopen(forecast_output_file, "r") do jld

            # read metadata
            metadata = read_forecast_metadata(jld)

            # read the DArray using read_darray
            data = DSGE.read_darray(jld)

            metadata, data
        end

        transforms, inds, date_inds = if pseudo_flag
            metadata[:pseudoobservable_revtransforms], metadata[:pseudoobservable_indices],
            metadata[:date_indices]
        elseif obs_flag
            metadata[:observable_revtransforms], metadata[:observable_indices], metadata[:date_indices]
        else
            error("means and bands are only calculated for observables and pseudoobservables")
        end

        # make sure date lists are valid
        date_list       = collect(keys(date_inds))   # array of actual dates
        date_inds_order = collect(values(date_inds)) # array of date indices
        check_consistent_order(date_list, date_inds_order)

        # make DataFrames for means and bands
        sort!(date_list)
        means = DataFrame(date = date_list)
        bands = Dict{Symbol,DataFrame}()

        # for each series (ie each pseudoobs, each obs, or each state):
        # 1. apply the appropriate transform
        # 2. add to DataFrame
        for (series, ind) in inds
            # apply transformation to all draws
            ex = Expr(:call, :map, parse_transform(transforms[series]), squeeze(data[:,ind,date_inds_order],2))
            show(dump(ex))
            transformed_data = eval(ex)

            # compute bands
            bands_one = find_density_bands(transformed_data, density_bands, minimize=false)
            bands_one[:date] = date_list

            # compute the mean and bands across draws and add to dataframe
            means[series] = vec(mean(transformed_data,1))
            bands[series] = bands_one
        end

        # put means and bands for this output_var in larger data structure
        if contains(string(output_var), "forecast")
            if pseudo_flag
                meanspseudo[:forecast] = means
                bandspseudo[:forecast] = bands
            elseif obs_flag
                meansobs[:forecast] = means
                bandsobs[:forecast] = bands
            end
        elseif contains(string(output_var), "hist")
            if pseudo_flag
                meanspseudo[:hist] = means
                bandspseudo[:hist] = bands
            elseif obs_flag
                meansobs[:hist] = means
                bandsobs[:hist] = bands
            end
        elseif contains(string(output_var), "shockdec")
            if pseudo_flag
                meanspseudo[:shockdec] = means
                bandspseudo[:shockdec] = bands
            elseif obs_flag
                meansobs[:shockdec] = means
                bandsobs[:shockdec] = bands
            end
        end
    end

    return MeansBands(:obs, meansobs, bandsobs), MeansBands(:pseudo, meanspseudo, bandspseudo)
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