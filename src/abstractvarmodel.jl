"""
```
AbstractVARModel{T} <: AbstractModel{T}
```

The AbstractVARModel is defined as a subtype of AbstractModel to implement
vector autoregression methods.
"""
abstract type AbstractVARModel{T} <: ModelConstructors.AbstractModel{T} end

function Base.show(io::IO, m::AbstractVARModel)
    @printf io "VAR Model\n"
    @printf io "observables:   %i\n" observables(m)
    @printf io "data_vintage:  %i\n" data_vintage(m)
    @printf io "description:\n %s\n" description(m)
end

get_parameters(m::AbstractVARModel) = m.parameters

"""
```
AbstractDSGEVARModel{T} <: AbstractVARModel{T}
```
The AbstractDSGEVARModel is defined as a subtype of AbstractVARModel to implement
DSGE-VAR methods.
"""
abstract type AbstractDSGEVARModel{T} <: AbstractVARModel{T} end

function Base.show(io::IO, m::AbstractDSGEVARModel)
    @printf io "DSGE-VAR Model\n"
    @printf io "observables:      %i\n" n_observables(m)
    @printf io "data_vintage:     %s\n" data_vintage(m)
    @printf io "DSGE model:       %s\n" spec(get_dsge(m))
    @printf io "DSGE description: %s\n" description(get_dsge(m))
end

"""
```
AbstractDSGEVECMModel{T} <: AbstractDSGEVARModel{T}
```
The AbstractDSGEVECMModel is defined as a subtype of AbstractDSGEVARModel to implement
DSGE-VECM methods.
"""
abstract type AbstractDSGEVECMModel{T} <: AbstractDSGEVARModel{T} end

# Number of anticipated policy shocks
n_anticipated_shocks(m::AbstractDSGEVARModel) = get_setting(m, :n_mon_anticipated_shocks)
n_anticipated_shocks_padding(m::AbstractDSGEVARModel) = get_setting(m, :n_mon_anticipated_shocks_padding)
n_mon_anticipated_shocks(m::AbstractDSGEVARModel) = get_setting(m, :n_mon_anticipated_shocks)
n_mon_anticipated_shocks_padding(m::AbstractDSGEVARModel) = get_setting(m, :n_mon_anticipated_shocks_padding)

# Dates, indices, number of periods for each regime
date_presample_start(m::AbstractDSGEVARModel) = get_setting(m, :date_presample_start)
date_mainsample_start(m::AbstractDSGEVARModel) = get_setting(m, :date_mainsample_start)
date_zlb_start(m::AbstractDSGEVARModel) = get_setting(m, :date_zlb_start)
# date_zlb_end(m::AbstractDSGEVARModel) = get_setting(m, :date_zlb_end)

date_presample_end(m::AbstractDSGEVARModel) = Dates.lastdayofquarter(get_setting(m, :date_mainsample_start) - Dates.Month(3))
date_prezlb_end(m::AbstractDSGEVARModel) = Dates.lastdayofquarter(get_setting(m, :date_zlb_start) - Dates.Month(3))
date_mainsample_end(m::AbstractDSGEVARModel) = Dates.lastdayofquarter(get_setting(m, :date_forecast_start) - Dates.Month(3))
date_conditional_end(m::AbstractDSGEVARModel) = get_setting(m, :date_conditional_end)

index_presample_start(m::AbstractDSGEVARModel) = 1
index_mainsample_start(m::AbstractDSGEVARModel) = subtract_quarters(date_mainsample_start(m), date_presample_start(m)) + 1
index_zlb_start(m::AbstractDSGEVARModel) = subtract_quarters(date_zlb_start(m), date_presample_start(m)) + 1
index_forecast_start(m::AbstractDSGEVARModel) = subtract_quarters(date_forecast_start(m), date_presample_start(m)) + 1

# Interface for fields and size
n_lags(m::AbstractDSGEVARModel)          = m.lags
get_lags(m::AbstractDSGEVARModel)        = m.lags
get_shocks(m::AbstractDSGEVARModel)      = m.shocks
get_λ(m::AbstractDSGEVARModel)           = m.λ
get_dsge(m::AbstractDSGEVARModel)        = m.dsge
n_observables(m::AbstractDSGEVARModel)   = length(m.observables)
n_cointegrating(m::AbstractDSGEVECMModel) = length(m.cointegrating)
n_cointegrating_add(m::AbstractDSGEVECMModel) = length(m.cointegrating_add)
n_shocks(m::AbstractDSGEVARModel)        = length(m.shocks)

# Interface for accessing parameters
get_parameters(m::AbstractDSGEVARModel) = m.dsge.parameters
n_parameters(m::AbstractDSGEVARModel)   = n_parameters(m.dsge)
n_parameters_steady_state(m::AbstractDSGEVARModel) = n_parameters_steady_state(m.dsge)
n_parameters_free(m::AbstractDSGEVARModel) = n_parameters_free(m.dsge)

# Interface for accessing rng
get_rng(m::AbstractDSGEVARModel) = m.dsge.rng

# Interface for accessing settings dictionary
get_settings(m::AbstractDSGEVARModel) = m.dsge.settings

# Interface for accessing observables dictionary
get_observables(m::AbstractDSGEVARModel) = m.observables

# Interface for accessing cointegrating dictionaries
get_cointegrating(m::AbstractDSGEVECMModel) = m.cointegrating
get_cointegrating_add(m::AbstractDSGEVECMModel) = m.cointegrating_add

# Interface for data
cond_vintage(m::AbstractDSGEVARModel)    = get_setting(m, :cond_vintage)
cond_id(m::AbstractDSGEVARModel)         = get_setting(m, :cond_id)
cond_full_names(m::AbstractDSGEVARModel) = get_setting(m, :cond_full_names)
cond_semi_names(m::AbstractDSGEVARModel) = get_setting(m, :cond_semi_names)
use_population_forecast(m::AbstractDSGEVARModel) = get_setting(m, :use_population_forecast)
hpfilter_population(m::AbstractDSGEVARModel)     = get_setting(m, :hpfilter_population)
n_conditional_periods(m::AbstractDSGEVARModel)   = n_conditional_periods(get_dsge(m))

# Interface for general computation settings
use_parallel_workers(m::AbstractDSGEVARModel)    = get_setting(m, :use_parallel_workers)

# Interface for estimation settings
reoptimize(m::AbstractDSGEVARModel)          = get_setting(m, :reoptimize)
calculate_hessian(m::AbstractDSGEVARModel) = get_setting(m, :calculate_hessian)
hessian_path(m::AbstractDSGEVARModel)      = get_setting(m, :hessian_path)
n_hessian_test_params(m::AbstractDSGEVARModel) = get_setting(m, :n_hessian_test_params)


# Interface for Metropolis-Hastings settings
n_mh_blocks(m::AbstractDSGEVARModel)       =  get_setting(m, :n_mh_blocks)
n_mh_param_blocks(m::AbstractDSGEVARModel) =  get_setting(m, :n_mh_param_blocks)
n_mh_simulations(m::AbstractDSGEVARModel)  =  get_setting(m, :n_mh_simulations)
n_mh_burn(m::AbstractDSGEVARModel)         =  get_setting(m, :n_mh_burn)
mh_thin(m::AbstractDSGEVARModel)           =  get_setting(m, :mh_thin)

# Interface for forecast settings
date_forecast_start(m::AbstractDSGEVARModel)   = get_setting(m, :date_forecast_start)
forecast_block_size(m::AbstractDSGEVARModel)   = get_setting(m, :forecast_block_size)
forecast_start_block(m::AbstractDSGEVARModel)  = get_setting(m, :forecast_start_block)
forecast_input_file_overrides(m::AbstractDSGEVARModel) = get_setting(m, :forecast_input_file_overrides)
forecast_uncertainty_override(m::AbstractDSGEVARModel) = get_setting(m, :forecast_uncertainty_override)
forecast_smoother(m::AbstractDSGEVARModel)     = get_setting(m, :forecast_smoother)
forecast_tdist_df_val(m::AbstractDSGEVARModel) = get_setting(m, :forecast_tdist_df_val)
forecast_tdist_shocks(m::AbstractDSGEVARModel) = get_setting(m, :forecast_tdist_shocks)
forecast_zlb_value(m::AbstractDSGEVARModel)    = get_setting(m, :forecast_zlb_value)
impulse_response_horizons(m::AbstractDSGEVARModel) = get_setting(m, :impulse_response_horizons)
n_shockdec_periods(m::AbstractDSGEVARModel)    = index_shockdec_end(m) - index_shockdec_start(m) + 1

# Interface for alternative policy settings
alternative_policy(m::AbstractDSGEVARModel) = alternative_policy(get_dsge(m))

function date_forecast_end(m::AbstractDSGEVARModel)
    date = date_forecast_start(m) + Dates.Month(3 * (forecast_horizons(m)-1))
    return Dates.lastdayofquarter(date)
end

function forecast_horizons(m::AbstractDSGEVARModel; cond_type::Symbol = :none)
    horizons = get_setting(m, :forecast_horizons)
    if cond_type == :none
        return horizons
    else
        return horizons - n_conditional_periods(m)
    end
end

function date_shockdec_start(m::AbstractDSGEVARModel)
    startdate = get_setting(m, :shockdec_startdate)
    if (Nullables.isnull(startdate)) | (startdate===nothing)
        return DSGEVAR.date_mainsample_start(m)
    else
        return get(startdate)
    end
end

function date_shockdec_end(m::AbstractDSGEVARModel)
    enddate =  get_setting(m, :shockdec_enddate)
    if !isnull(enddate)
        return get(enddate)
    else
        return date_forecast_end(m)
    end
end

"""
```
transform_to_model_space!(m::AbstractDSGEVARModel, values::Vector{T}) where T<:AbstractFloat
```

Transforms `values` from the real line to the model space, and assigns `values[i]` to
`m.parameters[i].value` for non-steady-state parameters. Recomputes the steady-state
paramter values.

### Arguments
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.
"""
function transform_to_model_space!(m::AbstractDSGEVARModel, values::Vector{T}) where {T<:AbstractFloat}
    new_values = transform_to_model_space(m.dsge.parameters, values)
    DSGE.update!(m, new_values)
    steadystate!(get_dsge(m))
end

# Updating parameters
function update!(m::AbstractDSGEVARModel{T}, values::Vector{T}) where {T<:Real}
    DSGE.update!(get_dsge(m), values)
    steadystate!(get_dsge(m))
end

function update!(m::AbstractDSGEVARModel{T}, values::ParameterVector{T}) where {T<:Real}
    DSGE.update!(get_dsge(m), values)
    steadystate!(get_dsge(m))
end

# Overload <= so that it applies to the underlying DSGE object
function (<=)(m::AbstractDSGEVARModel{T}, p::ModelConstructors.AbstractParameter{T}) where {T}
    dsge = get_dsge(m)
    dsge <= p
end

function (<=)(m::AbstractDSGEVARModel{T}, p::Union{ModelConstructors.SteadyStateParameter, ModelConstructors.SteadyStateParameterArray}) where {T}
    dsge = get_dsge(m)
    dsge <= p
end

function (<=)(m::AbstractDSGEVARModel{T}, p::ModelConstructors.SteadyStateParameterGrid) where {T}
    dsge = get_dsge(m)
    dsge <= p
end

function (<=)(m::AbstractDSGEVARModel, s::Setting)
    dsge = get_dsge(m)
    dsge <= s
end

# Setting access
function get_setting(m::AbstractDSGEVARModel, k::Symbol)
    return ModelConstructors.get_setting(get_dsge(m), k)
end

# Prior
function prior(m::AbstractDSGEVARModel)
    return ModelConstructors.prior(get_dsge(m))
end

# Not exposed to user. Actually create path and insert model string to file name.
function savepath(m::AbstractDSGEVARModel,
                  out_type::String,
                  sub_type::String,
                  file_name::String = "",
                  filestring_addl::Vector{String} = Vector{String}())
    # Containing directory
    dir = String(joinpath(saveroot(m), "output_data", spec(get_dsge(m)), subspec(get_dsge(m)),
                          spec(m), subspec(m), out_type, sub_type))

    if !isempty(file_name)
        base = filestring_base(m)
        return savepath(dir, file_name, base, filestring_addl)
    else
        return dir
    end
end

function savepath(dir::String,
                  file_name::String = "",
                  filestring_base::Vector{String} = Vector{String}(),
                  filestring_addl::Vector{String} = Vector{String}())
    if !isdir(dir)
        mkpath(dir)
    end

    if !isempty(file_name)
        (base, ext) = splitext(file_name)
        myfilestring = filestring(filestring_base, filestring_addl)
        file_name_detail = base * myfilestring * ext

        return joinpath(dir, file_name_detail)
    else
        return dir
    end
end

function filestring_base(m::AbstractDSGEVARModel)
    if !m.testing
        base = Vector{String}()
        for (skey, sval) in get_settings(m)
            if sval.print
                push!(base, ModelConstructors.to_filestring(sval))
            end
        end
        return base
    else
        return ["test"]
    end
end


# Helper to check for valid VAR systems
function check_valid_dsgevar(m::AbstractDSGEVARModel;
                             observables::Vector{Symbol} = collect(keys(m.observables)),
                             shocks::Vector{Symbol} = collect(keys(m.shocks)),
                             lags::Int = n_lags(m), λ::T = m.λ) where {T<:Real}
    # Check lags
    @assert lags >= 0 "Number of lags must be non-negative."

    # Check λ is non-negative
    @assert λ >= 0 "The weight on the DSGE prior λ must be non-negative."

    # Check observables
    missing_obs = Vector{Symbol}(undef, 0)
    for k in observables
        if !(k in keys(get_observables(get_dsge(m)))) &&
            !(k in keys(get_pseudo_observables(get_dsge(m))))
            push!(missing_obs, k)
        end
    end

    if !isempty(missing_obs)
        error("The following observables are not found in the underlying DSGE: \n" *
              join(string.(missing_obs), ", "))
    end

    # Check shocks
    missing_shocks = Vector{Symbol}(undef, 0)
    for k in shocks
        if !(k in keys(get_exogenous_shocks(get_dsge(m))))
            push!(missing_shocks, k)
        end
    end

    if !isempty(missing_shocks)
        error("The following shocks are not found in the underlying DSGE: \n" *
              join(string.(missing_shocks), ", "))
    end
end

# Helper to check for valid VECM systems
function check_valid_dsgevecm(m::AbstractDSGEVECMModel;
                              observables::Vector{Symbol} = collect(keys(m.observables)),
                              cointegrating::Vector{Symbol} = Vector{Symbol}(undef, 0),
                              cointegrating_add::Vector{Symbol} = Vector{Symbol}(undef, 0),
                              shocks::Vector{Symbol} = collect(keys(m.shocks)),
                              lags::Int = n_lags(m), λ::T = m.λ) where {T<:Real}

    # Check non-cointegrating fields
    check_valid_dsgevar(m; observables = observables,
                        shocks = shocks, lags = lags, λ = λ)


    # Check cointegrated observables
    missing_coint = Vector{Symbol}(undef, 0)
    for k in cointegrating
        if !(k in keys(get_observables(get_dsge(m)))) &&
            !(k in keys(get_pseudo_observables(get_dsge(m))))
            push!(missing_coint, k)
        end
    end

    if !isempty(missing_coint)
        error("The following cointegrating relationships are not found in the underlying DSGE: \n" *
              join(string.(missing_coint), ", "))
    end

    # Check additional cointegrated observables
    missing_coint = Vector{Symbol}(undef, 0)
    for k in cointegrating
        if !(k in keys(get_observables(get_dsge(m)))) &&
            !(k in keys(get_pseudo_observables(get_dsge(m))))
            push!(missing_coint, k)
        end
    end

    if !isempty(missing_coint)
        error("The following cointegrating relationships are not found in the underlying DSGE: \n" *
              join(string.(missing_coint), ", "))
    end
end

# Update VAR system, e.g. observables
function update!(m::AbstractDSGEVARModel; observables::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 shocks::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 lags::Int = 0, λ::T = m.λ, check_valid::Bool = true) where {T<:Real}

    if check_valid
        check_obs = isempty(observables) ? collect(keys(get_observables(m))) : observables
        check_sh  = isempty(shocks) ? collect(keys(get_shocks(m))) : shocks
        check_valid_dsgevar(m; observables = check_obs, shocks = check_sh,
                            lags = lags, λ = λ)
    end

    if m.λ != λ
        m.λ = λ
    end
    if !isempty(observables)
        for k in keys(m.observables)
            delete!(m.observables, k)
        end
        for (i,k) in enumerate(observables)
            m.observables[k] = i
        end
    end
    if !isempty(shocks)
        for k in keys(m.shocks)
            delete!(m.shocks, k)
        end
        for (i,k) in enumerate(shocks)
            m.shocks[k] = i
        end
    end
    if lags > 0
        m.lags = lags
    end

    return m
end

function update!(m::AbstractDSGEVARModel, dsge::AbstractDSGEModel;
                 observables::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 shocks::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 lags::Int = 0, λ::T = m.λ, check_valid::Bool = true) where {T<:Real}

    m.dsge = dsge
    update!(m; observables = observables, shocks = shocks, lags = lags,
            λ = λ, check_valid = check_valid)

    return m
end

# Update VECM system, e.g. observables
function update!(m::AbstractDSGEVECMModel; observables::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 cointegrating::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 cointegrating_add::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 shocks::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 lags::Int = 0, λ::T = m.λ, check_valid::Bool = true) where {T<:Real}

    if check_valid
        check_obs = isempty(observables) ? collect(keys(get_observables(m))) : observables
        check_sh  = isempty(shocks) ? collect(keys(get_shocks(m))) : shocks
        check_co  = isempty(cointegrating) ? collect(keys(get_cointegrating(m))) : cointegrating
        check_coa = isempty(cointegrating_add) ? collect(keys(get_cointegrating_add(m))) : cointegrating_add

        check_valid_dsgevecm(m; observables = check_obs,
                             cointegrating = check_co,
                             cointegrating_add = check_coa,
                             shocks = check_sh, lags = lags, λ = λ)
    end

    if m.λ != λ
        m.λ = λ
    end
    if !isempty(observables)
        for k in keys(get_observables(m))
            delete!(get_observables(m), k)
        end
        for (i,k) in enumerate(observables)
            m.observables[k] = i
        end
    end
    if !isempty(cointegrating)
        n_obs = length(m.observables)
        for k in keys(get_cointegrating(m))
            delete!(get_cointegrating(m), k)
        end
        for (i,k) in enumerate(cointegrating)
            m.cointegrating[k] = i + n_obs
        end
    end
    if !isempty(cointegrating_add)
        for k in keys(get_cointegrating_add(m))
            delete!(get_cointegrating_add(m), k)
        end
        for (i,k) in enumerate(cointegrating_add)
            m.cointegrating_add[k] = i
        end
    end
    if !isempty(shocks)
        for k in keys(m.shocks)
            delete!(m.shocks, k)
        end
        for (i,k) in enumerate(shocks)
            m.shocks[k] = i
        end
    end
    if lags > 0
        m.lags = lags
    end

    return m
end

function update!(m::AbstractDSGEVECMModel, dsge::AbstractDSGEModel;
                 observables::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 cointegrating::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 cointegrating_add::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 shocks::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 lags::Int = 0, λ::T = m.λ, check_valid::Bool = true) where {T<:Real}

    m.dsge = dsge
    update!(m; observables = observables, cointegrating = cointegrating,
            cointegrating_add = cointegrating_add,
            shocks = shocks, lags = lags,
            λ = λ, check_valid = check_valid)

    return m
end
