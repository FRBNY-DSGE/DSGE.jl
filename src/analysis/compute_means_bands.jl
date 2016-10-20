"""
```
type MeansBands
```

Stores the means and bands of results for a particular set of outputs from the forecast step.

Specifically, forecasts can be made for any element in the Cartesian product of 4 sets:

1. `input_type`: some subset of the parameter draws from the estimation step. See `forecast_all` for all possible options.
2. `cond_type`: conditional type. See `forecast_all` for all possible options.
3. *product*: a particular result computed in the forecast. This could be one of the following:
  - `hist`: smoothed histories
  - `forecast`: forecasted values
  - `shockdec`: shock decompositions
4. variable *class*: the category in which a particular variable, like `:y_t`, falls. Options are:
  - `state`: state (from `m.endogenous_states` or `m.endogenous_states_augmented`)
  - `obs`: observable (from `m.observables`)
  - `pseudo`: pseudoobservable (from `pseudo_measurement` equation)
  - `shock`: shock (from `m.exogenous_shocks`)

Note that the Cartesian product (product x class) is the set of options for `output_vars` in the `forecast_one` function signature.

### Fields

- `metadata::Dict{Symbol,Any}`: Contains metadata keeping track of the
  `input_type`, `cond_type`, product (history, forecast, shockdec,
  etc), and variable class (observable, pseudoobservable, state, etc)
  stored in this `MeansBands` structure.
- `means::DataFrame`: a `DataFrame` of the mean of the time series
- `bands::Dict{Symbol,DataFrame}`: a `Dict` mapping variable names to
  `DataFrame`s containing confidence bands for each variable. See
  `find_density_bands` for more information.
"""
type MeansBands
    metadata::Dict{Symbol,Any}
    means::DataFrame
    bands::Dict{Symbol,DataFrame}

    function MeansBands(key, means, bands)

        if !isempty(bands)
            # assert that means and bands fields have the same keys (provide info for same products)
            @assert sort(setdiff(names(means),[:date])) == sort(collect(keys(bands)))

            # check to make sure that # of periods in all dataframes are the same
            n_periods_means = size(means,1)
            for df in values(bands)
                n_periods_bands = size(df,1)
                @assert(n_periods_means == n_periods_bands,
                        "means and bands must have same number of periods")
            end
        end

        new(key, means, bands)
    end
end

"""
```
compute_means_bands{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])

compute_means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol, output_vars::Vector{Symbol},
                                               cond_type::Symbol, forecast_files::Dict{Symbol,S};
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])
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
"""
function compute_means_bands_all{T<:AbstractFloat}(m::AbstractModel, input_type::Symbol,
                                               output_vars::Vector{Symbol}, cond_type::Symbol;
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])

    # Get names of files that the forecast wrote
    forecast_output_files =  DSGE.get_output_files(m, input_type, output_vars, cond_type)

    compute_means_bands(input_type, output_vars, cond_type, forecast_output_files, density_bands = density_bands)
end
function compute_means_bands_all{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                               output_vars::Vector{Symbol},
                                               cond_type::Symbol,
                                               forecast_output_files::Dict{Symbol,S};
                                               density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])

    # set mb_output_vars so mb is the prefix to all meansbands output files
    mb_output_vars = [symbol("mb$x") for x in output_vars]

    mb_files = Dict{Symbol,UTF8String}()
    for x in keys(forecast_output_files)
        fn = forecast_output_files[x]
        dir  = dirname(fn)
        base = "mb" * basename(fn)
        mb_files[x] = joinpath(dir,base)
    end

    for (i,output_var) in enumerate(output_vars)

        # compute means and bands object
        mb = compute_means_bands(input_type, output_var, cond_type,
                                 forecast_output_files, density_bands = density_bands)

        # write to file
        filepath = mb_files[output_vars[i]]
        jldopen(filepath, "w") do file
            write(file, "mb", mb)
        end
        info("Wrote means and bands for $output_var to $filepath.")
    end
end

function compute_means_bands{T<:AbstractFloat, S<:AbstractString}(input_type::Symbol,
                                                    output_var::Symbol,
                                                    cond_type::Symbol,
                                                    forecast_output_files::Dict{Symbol,S};
                                                    density_bands::Array{T} = [0.5, 0.6, 0.7, 0.8, 0.9])

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

    # open correct input file
    forecast_output_file = forecast_output_files[output_var]
    metadata, data = jldopen(forecast_output_file, "r") do jld

        # read metadata
        metadata = read_forecast_metadata(jld)

        # read the DArray using read_darray
        data = DSGE.read_darray(jld)

        metadata, data
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
        transformed_data = eval(ex)

        # compute bands
        bands_one = find_density_bands(transformed_data, density_bands, minimize=false)
        bands_one[:date] = date_list

        # compute the mean and bands across draws and add to dataframe
        means[series] = vec(mean(transformed_data,1))
        bands[series] = bands_one
    end

    mb_metadata = Dict{Symbol,Any}(
                   :input_type => input_type,
                   :cond_type  => cond_type,
                   :product    => product,
                   :class      => class)

    return MeansBands(mb_metadata, means, bands)
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