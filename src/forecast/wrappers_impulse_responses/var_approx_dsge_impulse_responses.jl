"""
```
function impulse_responses(m::AbstractDSGEModel, paras::Union{Vector{S}, Matrix{S}},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_shock::Int; parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           use_intercept::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
```
computes the impulse responses of a VAR(p) approximation to a DSGE.

### Inputs
* `m::Union{AbstractDSGEModel,AbstractDSGEVARModel}`: DSGE/DSGEVAR model object
* `paras::Matrix{S}` or `paras::Vector{S}`: parameters to calibrate the model
* `input_type::Symbol`: `:mode` specifies a modal impulse response, and
    `:full` specifies a full-distribution forecast if `paras` is not given.
    This argument is also used to construct the file names of computed `MeansBands`.
* `method::Symbol`: type of impulse response to compute. The options are
    `:cholesky`, `:maximum_business_cycle_variance` or `:maxBC`,
    and `:cholesky_long_run` or `:choleskyLR`. See `?cholesky_shock`,
    `?maxBC_shock`, and `?choleskyLR_shock`.
* `lags::Int`: number of lags in the VAR(p) approximation, i.e. p = lags
* `observables::Vector{Symbol}`: observables to be used in the VAR. These can be
    any of the observables or pseudo-observables in `m`.
* `shocks::Vector{Symbol}`: (structural) exogenous shocks to be used in the DSGE-VAR.
    These shocks must be in `m`.
* `n_obs_shock::Int`: the index of the observable corresponding to the orthogonalized shock causing the impulse response.

### Keywords
* `parallel::Bool`: use parallel workers or not
* `frequency_band::Tuple{S,S}`: See `?maxBC_shock`.
* `flip_shocks::Bool`: impulse response shocks are negative by default. Set to `true` for
    a positive signed shock.
* `density_bands::Vector{Float64}`: bands for full-distribution IRF computations
* `create_meansbands::Bool`: set to `true` to save output as a `MeansBands` object.
* `minimize::Bool`: choose shortest interval if true, otherwise just chop off lowest and
    highest (percent/2)
* `forecast_string::String`: string tag for identifying this impulse response
* `verbose::Symbol`: quantity of output desired

"""
function impulse_responses(m::AbstractDSGEModel, paras::Union{Vector{S}, Matrix{S}},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_shock::Int; parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           use_intercept::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}

    # Compute VAR coefficients implied by DSGE
    if isa(m, AbstractDSGEModel)
        _dsgevar = DSGEVAR(m, shocks, "ss0")
        DSGE.update!(_dsgevar, lags = lags, observables = observables)
    else
        _dsgevar = m

    end

    return impulse_responses(_dsgevar, paras, input_type, method, n_obs_shock;
                             parallel = parallel, frequency_band = frequency_band,
                             flip_shocks = flip_shocks, use_intercept = use_intercept,
                             density_bands = density_bands, save_as_DSGE = isa(m, AbstractDSGEModel),
                             create_meansbands = create_meansbands, minimize = minimize,
                             forecast_string = forecast_string, verbose = verbose)
end

function impulse_responses(m::AbstractDSGEVARModel, paras::Union{Vector{S}, Matrix{S}},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_shock::Int; parallel::Bool = false,
                           use_intercept::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           save_as_DSGE::Bool = false,
                           verbose::Symbol = :high) where {S<:Real}
    DSGE.update!(m, lags = lags, observables = observables, shocks = shocks)
    return impulse_responses(m, paras, input_type, method, n_obs_shock;
                             parallel = parallel,
                             use_intercept = use_intercept,
                             frequency_band = frequency_band,
                             flip_shocks = flip_shocks,
                             density_bands = density_bands,
                             create_meansbands = create_meansbands,
                             minimize = minimize,
                             forecast_string = forecast_string,
                             save_as_DSGE = save_as_DSGE, verbose = verbose)
end


function impulse_responses(m::AbstractDSGEVARModel, paras::Union{Vector{S}, Matrix{S}},
                           input_type::Symbol, method::Symbol,
                           n_obs_shock::Int; parallel::Bool = false,
                           use_intercept::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false, test_meansbands::Bool = false,
                           minimize::Bool = true, save_as_DSGE::Bool = false,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    if n_obs_shock <= 0
        error("To use method $method, user must specify the index of" *
              " the target observable with keyword `n_obs_shock`.")
    end
    if isa(paras, Vector{S})
        paras = reshape(paras, 1, length(paras))
    end


    # Set up computation method
    mapfcn = parallel ? pmap : map
    h = impulse_response_horizons(m)

    # Compute IRFs
    paras = mapslices(x -> [vec(x)], paras, dims = 2)
    function get_irf!(para)
        DSGE.update!(m, para)
        return impulse_responses(m, method, n_obs_shock;
                                 horizon = h, use_intercept = use_intercept,
                                 flip_shocks = flip_shocks,
                                 frequency_band = frequency_band)
    end

    irf_output = mapfcn(para -> get_irf!(para), paras)

    if create_meansbands
        # Set up metadata and output from IRFs computation
        metadata = Dict{Symbol,Any}()
        metadata[:para] = input_type
        metadata[:cond_type] = :none
        metadata[:product] = :dsgevarirf
        metadata[:class] = :obs # We default everything to an observable
        metadata[:date_inds] = OrderedDict()

        # Set up for loop over variable names
        means = DataFrame()
        bands = Dict{Symbol,DataFrame}()
        observables = collect(keys(get_observables(m)))
        metadata[:indices] =
            OrderedDict{Symbol,Int}(name => name_i
                                    for (name_i, name) in enumerate(observables))

        # Means and Bands for each variable in a class
        for (name_i,name) in enumerate(observables)
            # irf_output is Vector{nobs x nperiod} -> for each observable,
            # we want to select its specific IRF, i.e. map(x -> x[obs_index, :]).
            # This creates a nperiod x ndraws matrix, which we want to transpose
            # to get a ndraws x nperiod matrix
            single_var = Matrix(reduce(hcat, map(x -> x[name_i, :], irf_output))')
            @show size(single_var)
            means[!,name] = vec(mean(single_var, dims = 1))
            bands[name]   = find_density_bands(single_var, density_bands;
                                               minimize = minimize)
        end
        mb = MeansBands(metadata, means, bands)

        # Save MeansBands
        if !test_meansbands
            tail = if method == :cholesky
                :cholesky
            elseif method == :maxBC || method == :maximum_business_cycle_variance
                :maxBC
            else
                :choleskyLR
            end

            var_names = save_as_DSGE ?
                Symbol("_" * join(string.(map(x -> DSGE.detexify(x), observables)), "_") * "_") : Symbol("_")
            fp = get_meansbands_output_file(save_as_DSGE ? get_dsge(m) : m, input_type, :none,
                                            Symbol(:dsgevarirf, :obs,
                                                   var_names, tail),
                                            forecast_string = forecast_string)
            dirpath = dirname(fp)
            isdir(dirpath) || mkpath(dirpath)
            JLD2.jldopen(fp, true, true, true, IOStream) do file
                write(file, "mb", mb)
            end
            println(verbose, :high, "  " * "wrote " * basename(fp))
        end
        return mb
    else
        # Reshape irf_output to nobs x nperiod x ndraw
        # return cat(map(x -> x', irf_output)..., dims = 3)
        return cat(irf_output..., dims = 3)
    end
end
