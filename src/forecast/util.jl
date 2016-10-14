"""
```
compute_system(m)
```

Given the current model parameters, compute the state-space system
corresponding to model `m`. Returns a `System` object.
"""
function compute_system{T<:AbstractFloat}(m::AbstractModel{T})
    # Solve model
    TTT, RRR, CCC = solve(m)
    transition_equation = Transition(TTT, RRR, CCC)

    # Solve measurement equation
    shocks = n_anticipated_shocks(m) > 0
    measurement_equation = measurement(m, TTT, RRR, CCC; shocks = shocks)

    # Solve pseudo-measurement equation
    pseudo_measurement_equation = if method_exists(pseudo_measurement, (typeof(m),)) && forecast_pseudoobservables(m)
        _, pseudo_mapping = pseudo_measurement(m)
        Nullable(pseudo_mapping)
    else
        Nullable{PseudoObservableMapping{T}}()
    end

    return System(transition_equation, measurement_equation, pseudo_measurement_equation)
end

"""
```
get_jstep(m, n_sim)
```

Retrieve `forecast_jstep` setting (thinning step size for forecast
step) from `m.settings`. If `n_sim ==  1`, returns 1.
"""
function get_jstep(m, n_sim)
    if n_sim == 1
        jstep = 1
    else
        jstep = get_setting(m, :forecast_jstep)
    end
end

"""
```
get_input_file(m, input_type)
```

Compute the appropriate forecast input filenames for model `m` and
forecast input type `input_type`.

The default input files for each `input_type` can be overriden by adding entries
to the `Dict{Symbol, ASCIIString}` returned from
`forecast_input_file_overrides(m)`. For example:

```
overrides = forecast_input_file_overrides(m)
overrides[:mode] = \"path/to/input/file.h5\"
```
"""
function get_input_file(m, input_type)
    overrides = forecast_input_file_overrides(m)
    if haskey(overrides, input_type)
        override_file = overrides[input_type]
        if ispath(override_file)
            return override_file
        else
            error("Invalid input file override for input_type = $input_type: $override_file")
        end
    end

    if input_type == :mode
        return rawpath(m,"estimate","paramsmode.h5")
    elseif input_type == :mean
        return workpath(m,"estimate","paramsmean.h5")
    elseif input_type == :init
        return ""
    elseif input_type == :full
        return rawpath(m,"estimate","mhsave.h5")
    elseif input_type == :subset
        throw(ArgumentError("Not implemented."))
    else
        throw(ArgumentError("Invalid input_type: $(input_type)"))
    end
end

"""
```
get_output_vars(m, output_type)
```

Returns a vector of `output_vars` corresponding to the `output_type`. See the
`forecast_all` documentation for a complete list of possible `output_type`s.
"""
function get_output_vars(m, output_type)
    if output_type == :states
        vars = [:histstates,
                :histpseudo]
    elseif output_type == :shocks
        vars = [:histshocks]
    elseif output_type == :shocks_nonstandardized
        vars = [:histshocksns]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :forecast
       vars = [:forecaststates,
               :forecastobs,
               :forecastpseudo,
               :forecastshocks]
    elseif output_type == :shockdec
        vars = [:shockdecstates,
                :shockdecpseudo,
                :shockdecobs]
    elseif output_type == :dettrend
        vars = [:dettrendstates,
                :dettrendobs]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :counter
        vars = [:counterstates,
                :counterobs]
        throw(ArgumentError("Not implemented."))
    elseif output_type == :simple
        vars = [:histstates,
                :histpseudo,
                :forecaststates,
                :forecastpseudo,
                :forecastobs,
                :forecastshocks]
    elseif output_type == :all
        vars = []
        throw(ArgumentError("Not implemented."))
    else
        throw(ArgumentError("Invalid input_type: $(output_type)"))
    end
end

"""
```
get_output_files(m, input_type, output_vars, cond_type)
```

Compute the appropriate forecast output filenames for model `m`, forecast input
type `input_type`, and conditional type `cond_type`, for each output variable in
`output_vars`.
"""
function get_output_files(m, input_type, output_vars, cond_type)
    additional_file_strings = ASCIIString[]
    push!(additional_file_strings, "para=" * abbrev_symbol(input_type))
    push!(additional_file_strings, "cond=" * abbrev_symbol(cond_type))

    return [symbol(x) => rawpath(m, "forecast", "$x.jld", additional_file_strings) for x in output_vars]
end

typealias DVector{T, A} DArray{T, 1, A}
typealias DMatrix{T, A} DArray{T, 2, A}

"""
```
write_darray(filepath, darr)
```

Write the contents of the `DArray` `darr` to the JLD file `filepath` without
converting `darr` back to an `Array` or copying local parts to the originator
process. `darr` can then be read back in as an `Array` using `read_darray`.

Each worker process `pid` writes its own local part and local indices to the
file as `arr\$pid` and `inds\$pid`. The dimensions `dims` of `darr` and the
vector of processes `pids` over which `darr` is distributed are also written to
`filepath` in order to facilitate reading back in.
"""
function write_darray{T<:AbstractFloat}(filepath::AbstractString, darr::DArray{T})
    function write_localpart(pid::Int)
        jldopen(filepath, "r+") do file
            write(file, "inds$pid", collect(localindexes(darr)))
            write(file, "arr$pid", localpart(darr))
        end
    end

    jldopen(filepath, "w") do file
        write(file, "dims", darr.dims)
        write(file, "pids", collect(darr.pids))
    end

    for pid in darr.pids
        remotecall_wait(pid, write_localpart, pid)
        sleep(0.001) # Need this, or else HDF5 complains that the new object already exists
    end
end

"""
```
read_darray(filepath)
```

Read the `DArray` saved in `filepath` by `write_darray`. Returns an `Array`.
"""
function read_darray(filepath::AbstractString)
    file = jldopen(filepath, "r")
    dims = read(file, "dims")
    pids = read(file, "pids")

    out = zeros(dims...)
    for pid in pids
        inds = read(file, "inds$pid")
        out[inds...] = read(file, "arr$pid")
    end
    close(file)
    return out
end

"""
```
write_forecast_metadata(m, filepath, var)
```

Write metadata about the saved forecast output `var` to `filepath`.

Specifically, we save dictionaries mapping dates, as well as state, observable,
pseudo-observable, and shock names, to their respective indices in the saved
forecast output array. The saved dictionaries include:

- `date_indices::Dict{Date, Int}`: saved for all forecast outputs
- `state_names::Dict{Symbol, Int}`: saved for `var in [:histstates, :forecaststates, :shockdecstates]`
- `observable_names::Dict{Symbol, Int}`: saved for `var in [:forecastobs, :shockdecobs]`
- `pseudoobservable_names::Dict{Symbol, Int}`: saved for `var in [:histpseudo, :forecastpseudo, :shockdecpseudo]`
- `shock_names::Dict{Symbol, Int}`: saved for `var in [:histshocks, :forecastshocks, :shockdecstates, :shockdecobs, :shockdecpseudo]`
"""
function write_forecast_metadata(m::AbstractModel, filepath::AbstractString, var::Symbol)
    file = jldopen(filepath, "r+")

    # Write date range
    dates = if contains(string(var), "hist")
        quarter_range(date_prezlb_start(m), date_zlb_end(m))
    elseif contains(string(var), "forecast")
        quarter_range(date_forecast_start(m), date_forecast_end(m))
    elseif contains(string(var), "shockdec")
        date_shockdec_start = if isnull(shockdec_startdate(m))
            date_prezlb_start(m)
        else
            get(shockdec_startdate(m))
        end
        date_shockdec_end = if isnull(shockdec_enddate(m))
            date_forecast_end(m)
        else
            get(shockdec_enddate(m))
        end
        quarter_range(date_shockdec_start, date_shockdec_end)
    end
    date_indices = [d::Date => i::Int for (i, d) in enumerate(dates)]
    write(file, "date_indices", date_indices)

    # Write state names
    if contains(string(var), "states")
        state_indices = merge(m.endogenous_states, m.endogenous_states_augmented)
        @assert length(state_indices) == n_states_augmented(m) # assert no duplicate keys
        write(file, "state_indices", state_indices)

    # Write observable names
    elseif contains(string(var), "obs")
        write(file, "observable_indices", m.observables)

    # Write pseudo-observable names
    elseif contains(string(var), "pseudo")
        _, pseudo_mapping = pseudo_measurement(m)
        write(file, "pseudoobservable_indices", pseudo_mapping.inds)

    # Write shock names
    elseif contains(string(var), "shocks") || contains(string(var), "shockdec")
        write(file, "shock_indices", m.exogenous_shocks)
    end

    close(file)
end