"""
`compute_system(m)`

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
`get_jstep(m, n_sim)`

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
`get_input_file(m, input_type)`

Computes the appropriate forecast input filenames for model `m` and
forecast input type `input_type`. For example, the default file containing the parameter mode is ``

Default input file
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


function get_output_files(m, input_type, output_vars, cond_type)
    additional_file_strings = ASCIIString[]
    push!(additional_file_strings, "para=" * abbrev_symbol(input_type))
    push!(additional_file_strings, "cond=" * abbrev_symbol(cond_type))

    return [symbol(x) => rawpath(m, "forecast", "$x.jld", additional_file_strings) for x in output_vars]
end


"""
`DVector{T, A}`

typalias for `DArray{T, 1, A}`
"""
typealias DVector{T, A} DArray{T, 1, A}

"""
`DMatrix{T, A}`

typealias for `DArray{T, 2, A}`
"""
typealias DMatrix{T, A} DArray{T, 2, A}


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
        sleep(0.001)
    end
end


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
