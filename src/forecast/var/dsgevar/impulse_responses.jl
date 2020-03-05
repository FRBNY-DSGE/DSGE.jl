"""
```
function impulse_responses(m, input_type, method,
                           lags, observables, shocks,
                           n_obs_var; parallel = false,
                           frequency_band = (2*π/32, 2*π/6),
                           flip_shocks = false,
                           density_bands = [.5, .6, .7, .8, .9],
                           create_meansbands = false,
                           minimize = true,
                           forecast_string = "",
                           verbose = :high) where {S<:Real}
function impulse_responses(m, paras, input_type, method,
                           lags, observables, shocks,
                           n_obs_var; parallel = false,
                           frequency_band = (2*π/32, 2*π/6),
                           flip_shocks = false,
                           density_bands = [.5, .6, .7, .8, .9],
                           create_meansbands = false,
                           minimize = true,
                           forecast_string = "",
                           verbose = :high) where {S<:Real}
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
* `n_obs_var::Int`: the index of the observable to be shocked by
    the reduced-form impulse response to the VAR system.

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

function impulse_responses(m::Union{AbstractDSGEModel,AbstractDSGEVARModel}, paras::Matrix{S},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_var::Int;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           use_intercept::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    if n_obs_var <= 0
        error("To use method $method, user must specify the index of" *
              " the target observable with keyword `n_obs_var`.")
    end

    # Set up computation method
    mapfcn = parallel ? pmap : map
    h = impulse_response_horizons(m)

    # Compute VAR coefficients implied by DSGE
    paras = mapslices(x -> [vec(x)], paras, dims = 2)
    get_β_Σ! = if isa(m, AbstractDSGEVARModel)
        function _get_β_Σ_dsgevar_(para)
            DSGE.update!(m, para)
            return compute_system(m; verbose = verbose, use_intercept = use_intercept)
        end
    else
        use_measurement_error = hasmethod(measurement_error, (typeof(m),))

        function _get_β_Σ_dsge_(para)
            DSGE.update!(m, para)
            system = compute_system(m; verbose = verbose)
            system = compute_system(dsge, system; observables = observables, shocks = shocks)
            nobs = length(observables)
            nshocks = length(shocks)
            EE, MM = use_measurement_error ? measurement_error(m) :
                zeros(nobs, nobs), zeros(nobs, nshocks)
            return var_approx_state_space(system[:TTT], system[:RRR], system[:QQ],
                                          system[:DD], system[:ZZ], EE, MM, lags;
                                          get_population_moments = false)
        end
    end

    var_output =
        mapfcn(para -> get_β_Σ!(para), paras)

    # Reformat output
    β_draws = map(x -> x[1], var_output)
    Σ_draws = map(x -> x[2], var_output)

    # Compute IRFs
    irf_output =
        mapfcn((β, Σ) ->
               impulse_responses(β, Σ, n_obs_var, h; method = method,
                                 use_intercept = false,
                                 flip_shocks = flip_shocks),
               β_draws, Σ_draws)

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
        metadata[:indices] =
            OrderedDict{Symbol,Int}(name => name_i
                                    for (name_i, name) in enumerate(observables))

        # Means and Bands for each variable in a class
        for (name_i,name) in enumerate(observables)
            # irf_output is Vector{nperiod x nobs} -> for each observable,
            # we want to select its specific IRF, i.e. map(x -> x[:,obs_index]).
            # This creates a nperiod x ndraws matrix, which we want to transpose
            # to get a ndraws x nperiod matrix
            single_var = Matrix(reduce(hcat, map(x -> x[:,name_i], irf_output))')
            means[!,name] = vec(mean(single_var, dims = 1))
            bands[name]   = find_density_bands(single_var, density_bands;
                                               minimize = minimize)
        end
        mb = MeansBands(metadata, means, bands)

        # Save MeansBands
        tail = if method == :cholesky
            :cholesky
        elseif method == :maxBC || method == :maximum_business_cycle_variance
            :maxBC
        else
            :choleskyLR
        end

        var_names = isa(m, AbstractDSGEModel) ?
            Symbol("_" * join(string.(DSGE.detexify(observables)), "_") * "_") : Symbol("_")
        fp = get_meansbands_output_file(m, input_type, :none,
                                        Symbol(:dsgevarirf, :obs,
                                               var_names, tail),
                                        forecast_string = forecast_string)
        dirpath = dirname(fp)
        isdir(dirpath) || mkpath(dirpath)
        JLD2.jldopen(fp, true, true, true, IOStream) do file
            write(file, "mb", mb)
        end
        println(verbose, :high, "  " * "wrote " * basename(fp))
        return mb
    else
        # Reshape irf_output to nobs x nperiod x ndraw
        return cat(map(x -> x', irf_output)..., dims = 3)
    end
end

function impulse_responses(m::AbstractDSGEVARModel, paras::Vector{S},
                           input_type::Symbol, method::Symbol,
                           n_obs_var::Int; parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    return impulse_responses(m, reshape(paras, 1, length(paras)),
                             input_type, method, get_lags(m), collect(keys(get_observables(m))),
                             collect(keys(get_shocks(m))), n_obs_var;
                             parallel = parallel,
                             frequency_band = frequency_band,
                             flip_shocks = flip_shocks,
                             density_bands = density_bands,
                             create_meansbands = create_meansbands,
                             minimize = minimize,
                             forecast_string = forecast_string,
                             verbose = verbose)
end

function impulse_responses(m::AbstractDSGEModel, paras::Vector{S},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_var::Int; parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    return impulse_responses(m, reshape(paras, 1, length(paras)),
                             input_type, method, lags, observables,
                             shocks, n_obs_var;
                             parallel = parallel,
                             frequency_band = frequency_band,
                             flip_shocks = flip_shocks,
                             density_bands = density_bands,
                             create_meansbands = create_meansbands,
                             minimize = minimize,
                             forecast_string = forecast_string,
                             verbose = verbose)
end

function impulse_responses(m::AbstractDSGEVARModel, paras::Matrix{S},
                           input_type::Symbol, method::Symbol,
                           n_obs_var::Int; parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    return impulse_responses(m, paras, input_type, method,
                             get_lags(m), collect(keys(get_observables(m))),
                             collect(keys(get_shocks(m))), n_obs_var;
                             parallel = parallel,
                             frequency_band = frequency_band,
                             flip_shocks = flip_shocks,
                             density_bands = density_bands,
                             create_meansbands = create_meansbands,
                             minimize = minimize,
                             forecast_string = forecast_string,
                             verbose = verbose)
end

"""
```
function impulse_responses(m, input_type, method,
                           lags, observables, shocks,
                           n_obs_var; parallel = false,
                           frequency_band = (2*π/32, 2*π/6),
                           flip_shocks = false,
                           density_bands = [.5, .6, .7, .8, .9],
                           create_meansbands = false,
                           minimize = true,
                           forecast_string = "",
                           verbose = :high) where {S<:Real}
function impulse_responses(m, paras, input_type, method,
                           lags, observables, shocks,
                           n_obs_var; parallel = false,
                           frequency_band = (2*π/32, 2*π/6),
                           flip_shocks = false,
                           density_bands = [.5, .6, .7, .8, .9],
                           create_meansbands = false,
                           minimize = true,
                           forecast_string = "",
                           verbose = :high) where {S<:Real}
```
computes the VAR impulse responses identified by the state space system.

****
NOTE: this function generally involves taking random draws from
probability distributions, so seeds need to be set
to achieve reproducibility.
****

```
sₜ = TTT × sₜ₋₁ + RRR × impact[:, i],
yₜ = ZZ × sₜ + DD + MM × impact[:, i],
```
where `impact[:, i]` is a linear combination of
(orthogonal) structural shocks `ϵₜ ∼ 𝒩 (0, I)`, and
`MM × impact[:, i]` are the correlated measurement errors.

The VAR impulse responses are computed according to
```
ŷₜ₊₁ = X̂ₜ₊₁β + uₜ₊₁,
```
where `X̂ₜ₊₁` are the lags of observables in period `t + 1`, i.e. `yₜ, yₜ₋₁, ..., yₜ₋ₚ`.

If the method is `:rotation`, the shock `uₜ₊₁` is identified via
```
Σᵤ = 𝔼[u × u'] = chol(Σᵤ) × Ω × ϵₜ,
```
where the rotation matrix `Ω` is the `Q` matrix from a QR decomposition
of the impact response matrix corresponding to the state space system, i.e.
```
Ω, _ = qr(∂yₜ / ∂ϵₜ').
```

Otherwise, we draw a β and Σᵤ from the posterior implied by the DSGE
and data, and we then compute normal VAR impulse responses given those
coefficients and error covariance matrix.

### Inputs
* `m::Union{AbstractDSGEModel,AbstractDSGEVARModel}`: DSGE/DSGEVAR model object
* `paras::Matrix{S}` or `paras::Vector{S}`: parameters to calibrate the model
* `input_type::Symbol`: `:mode` specifies a modal impulse response, and
    `:full` specifies a full-distribution forecast if `paras` is not given.
    This argument is also used to construct the file names of computed `MeansBands`.
* `method::Symbol`: type of impulse response to compute. The options are
    `:cholesky` and `:rotation`. For the first, see `?cholesky_shock`,
    and for the latter, we use the DSGE model to identify the rotation matrix
    which maps the DSGE's structural shocks to the innovations in the VAR's observables.
* `lags::Int`: number of lags in the VAR(p) approximation, i.e. p = lags
* `observables::Vector{Symbol}`: observables to be used in the VAR. These can be
    any of the observables or pseudo-observables in `m`.
* `shocks::Vector{Symbol}`: (structural) exogenous shocks to be used in the DSGE-VAR.
    These shocks must be in `m`.
* `n_obs_var::Int`: the index of the observable to be shocked by
    the reduced-form impulse response to the VAR system.

### Keywords
* `parallel::Bool`: use parallel workers or not
* `frequency_band::Tuple{S,S}`: See `?maxBC_shock`.
* `n_obs_var::Int`: Index of observable to be shocked when using a Cholesky-based impulse response
* `draw_shocks::Bool`: true if you want to draw shocks along the entire horizon
* `flip_shocks::Bool`: impulse response shocks are negative by default. Set to `true` for
    a positive signed shock.
* `accumulate::Bool`: accumulate any impulse responses, must set `cum_inds` as well
* `cum_inds::Union{Int,UnitRange{int},Vector{Int}}`: indices of impulse responses to accumulate
* `density_bands::Vector{Float64}`: bands for full-distribution IRF computations
* `create_meansbands::Bool`: set to `true` to save output as a `MeansBands` object.
* `minimize::Bool`: choose shortest interval if true, otherwise just chop off lowest and
    highest (percent/2)
* `forecast_string::String`: string tag for identifying this impulse response
* `verbose::Symbol`: quantity of output desired

"""
function impulse_responses(m::AbstractDSGEVARModel{S}, paras::Matrix{S},
                           data::Matrix{S}, input_type::Symbol, method::Symbol;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           n_obs_var::Int = 1, draw_shocks::Bool = false,
                           flip_shocks::Bool = false, accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    @assert !isempty(data) "A non-empty data matrix with dimensions (n_observables x n_periods) must be passed to use $(string(method))"

    if !(method in [:rotation, :cholesky, :cholesky_long_run, :choleskyLR, :maxBC,
                      :maximum_business_cycle_variance])
        error("A VAR IRF identified by a DSGE's rotation using method"
              * " $(string(method)) has not been implemented.")
    end

    # Compute dimensions of needed objects
    lags = get_lags(m)
    T = size(data, 2) - lags
    nobs = size(data, 1)
    λ = get_λ(m)
    T_λT = convert(Int, T + λ * T)
    k = nobs * lags + 1
    h = impulse_response_horizons(m)

    # Get population moments
    YYYY, XXYY, XXXX =
        compute_var_population_moments(data, lags; use_intercept = true)
    if method == :rotation
        XX = lag_data(data, lags; use_intercept = true)
        X̂ = vcat(1, data[:, end], XX[end, 1+1:k - nobs])
    end

    # Get measurement error
    if hasmethod(measurement_error, (typeof(m),))
        _, MM = measurement_error(m)
    else
        MM = zeros(S, n_observables(m), n_shocks(m))
    end

    dsgevarrotationirf_method = if method == :rotation
        function _dsgevar_λ_rotation_irf_(para)
            DSGE.update!(m, para)
            sys = compute_system(m, get_system = true, use_intercept = true)
            β, Σ = compute_system(m, data; verbose = verbose, use_intercept = true)
            Σ += Σ'
            Σ ./= 2.

            return DSGE.impulse_responses(sys[:TTT], sys[:RRR], sys[:ZZ], sys[:DD], MM,
                                          sys[:QQ], k, β, Σ, X̂, h;
                                          method = method,
                                          accumulate = accumulate,
                                          cum_inds = cum_inds, flip_shocks = flip_shocks,
                                          draw_shocks = draw_shocks)
        end
    elseif method in [:cholesky, :cholesky_long_run, :choleskyLR, :maxBC,
                      :maximum_business_cycle_variance]
        function _dsgevar_λ_irf_(para)
            DSGE.update!(m, para)

            β, Σ = compute_system(m, data; verbose = verbose, use_intercept = true)
            Σ += Σ'
            Σ ./= 2.

            return impulse_responses(β, Σ, n_obs_var, h; method = method,
                              use_intercept = true,
                              flip_shocks = flip_shocks)
        end
    end

    mapfcn = parallel ? pmap : map
    paras = mapslices(x -> [vec(x)], paras, dims = 2)
    irf_output =
        mapfcn(para -> dsgevarrotationirf_method(para), paras)

    if create_meansbands
        # Set up metadata and output from IRFs computation
        metadata = Dict{Symbol,Any}()
        metadata[:para] = input_type
        metadata[:cond_type] = :none
        metadata[:product] = :dsgevarlambdairf
        metadata[:class] = :obs # We default everything to an observable
        metadata[:date_inds] = OrderedDict()

        # Set up for loop over variable names
        means = DataFrame()
        bands = Dict{Symbol,DataFrame}()
        metadata[:indices] = get_observables(m)

        # Means and Bands for each variable in a class
        for (name, name_i) in get_observables(m)
            # irf_output is Vector{nperiod x nobs} -> for each observable,
            # we want to select its specific IRF, i.e. map(x -> x[:,obs_index]).
            # This creates a nperiod x ndraws matrix, which we want to transpose
            # to get a ndraws x nperiod matrix
            single_var = Matrix(reduce(hcat, map(x -> x[:,name_i], irf_output))')
            means[!,name] = vec(mean(single_var, dims = 1))
            bands[name]   = find_density_bands(single_var, density_bands;
                                               minimize = minimize)
        end
        mb = MeansBands(metadata, means, bands)

        # Save MeansBands
        tail = if method == :cholesky
            :cholesky
        elseif method == :maxBC || method == :maximum_business_cycle_variance
            :maxBC
        elseif method == :rotation
            :rotation
        else
            :choleskyLR
        end

        fp = get_meansbands_output_file(m, input_type, :none,
                                        Symbol(metadata[:product], :obs_, tail),
                                        forecast_string = forecast_string)
        dirpath = dirname(fp)
        isdir(dirpath) || mkpath(dirpath)
        JLD2.jldopen(fp, true, true, true, IOStream) do file
            write(file, "mb", mb)
        end
        println(verbose, :high, "  " * "wrote " * basename(fp))
        return mb
    else
        # Reshape irf_output to nobs x nperiod x ndraw
        return cat(map(x -> x', irf_output)..., dims = 3)
    end
end

function impulse_responses(m::AbstractDSGEVARModel{S}, paras::Vector{S},
                           data::Matrix{S}, input_type::Symbol, method::Symbol;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           n_obs_var::Int = 1, draw_shocks::Bool = false,
                           flip_shocks::Bool = false, accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    return impulse_responses(m, reshape(paras, 1, length(paras)),
                             data, input_type, method; parallel = parallel,
                             frequency_band = frequency_band, n_obs_var = n_obs_var,
                             draw_shocks = draw_shocks, flip_shocks = flip_shocks,
                             accumulate = accumulate, cum_inds = cum_inds,
                             density_bands = density_bands,
                             create_meansbands = create_meansbands, minimize = minimize,
                             forecast_string = forecast_string, verbose = verbose)
end

"""
```
function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, impact::Matrix{S}, horizon::Int;
                           accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0) where {S<:Real}
```
computes the impulse responses of the linear state space system to linear combinations
of (orthogonal) structural shocks specified by `impact`.
Measurement error that is correlated with
the impact matrix is allowed. We also include the option to accumulate certain observables.

This state space model takes the form

```
sₜ = TTT × sₜ₋₁ + RRR × impact[:, i],
yₜ = ZZ × sₜ + DD + MM × impact[:, i],
```
where `impact[:, i]` is a linear combination of orthogonal structural shocks
with mean zero and identity covariance, and `MM × impact[:, i]` are the
correlated measurement errors.

The `impact` matrix must be `nshock × nirf`, where `nshock` is
 the number of structural shocks and `nirf` is the number
of desired impulse responses. For each row of `impact`,
we compute the corresponding impulse responses.

A standard choice for `impact` is a square diagonal matrix. In this case,
we compute the impulse response of observables to each structural shock,
scaled by the desired size of the shock.

### Keywords
* `accumulate`: set to true if an observable should be accumulated.
* `cum_inds`: specifies the indices of variables to accumulated.

### Outputs
* `irf_results::Matrix{S}`: a `nobs x horizon × nirf` matrix, where
     `nobs` is the number of observables.
"""
function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, impact::Matrix{S}, horizon::Int;
                           accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0) where {S<:Real}
    # Get dimensions
    nshock, nirf = size(impact)
    nstate       = size(TTT, 1)
    nobs         = size(ZZ, 1)

    # Compute impulse response to identified impact matrix
    irf_results = Array{S, 3}(undef, nobs, horizon, nirf)
    for i = 1:nirf
        imp    = impact[:, i]
        states = zeros(S, nstate, horizon)
        obs    = zeros(S, nobs, horizon)

        states[:, 1] = RRR * imp
        obs[:, 1]    = ZZ * states[:, 1] + MM * imp
        for t = 2:horizon
            states[:, t] = TTT * states[:, t - 1]
            obs[:, t]    = ZZ * states[:, t] + DD
        end
        if accumulate
            obs[cum_inds, :] = cumsum(obs[cum_inds, :], dims = length(cum_inds) > 1 ? 2 : 1)
        end
        irf_results[:, :, i] = obs
    end

    return irf_results
end

"""
```
function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, QQ::Matrix{S},
                           k::Int, β::Matrix{S}, Σ::Matrix{S},
                           x̂::Matrix{S}, horizon::Int;
                           accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0,
                           test_shocks::Matrix{S} =
                           Matrix{S}(undef, 0, 0)) where {S<:Real}
```
computes the VAR impulse responses identified by the state space system
```
sₜ = TTT × sₜ₋₁ + RRR × ϵₜ
yₜ = ZZ × sₜ + DD + MM × ϵₜ
```
where `ϵₜ ∼ 𝒩 (0, QQ)` and `MM × ϵₜ` are the correlated measurement errors.

The VAR impulse responses are computed according to
```
ŷₜ₊₁ = X̂ₜ₊₁β + uₜ₊₁,
```
where `X̂ₜ₊₁` are the lags of observables in period `t + 1`, i.e. `yₜ, yₜ₋₁, ..., yₜ₋ₚ`.
The shock `uₜ₊₁` is identified via
```
Σᵤ = 𝔼[u × u'] = chol(Σᵤ) × Ω × ϵₜ,
```
where the rotation matrix `Ω` is the `Q` matrix from a QR decomposition
of the impact response matrix corresponding to the state space system, i.e.
```
Ω, _ = qr(∂yₜ / ∂ϵₜ').
```
"""

function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, QQ::Matrix{S},
                           k::Int, β::Matrix{S}, Σ::Matrix{S},
                           X̂::Vector{S}, horizon::Int; method::Symbol = :rotation,
                           accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0,
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
                           test_shocks::Matrix{S} =
                           Matrix{S}(undef, 0, 0)) where {S<:Real}
    if !(method in [:rotation])
        error("Impulse responses for method $(string(method)) have not been implemented.")
    end

    # Compute impulse responses of predicted values for each β, Σ, and rotation
    nobs = size(ZZ, 1)
    a0_m = convert(Matrix{S},
                   dropdims(impulse_responses(TTT, RRR, ZZ, DD, MM,
                                              sqrt.(QQ), 1, # 1 b/c just want impact and sqrt -> 1 sd shock
                                              accumulate = accumulate,
                                              cum_inds = cum_inds); dims = 2)')
    rotation, _ = qr(a0_m)
    Σ_chol = cholesky(Σ).L * rotation' # mapping from structural shocks to innovations in VAR

    if draw_shocks || !isempty(test_shocks)
        ŷ = Matrix{S}(undef, horizon, nobs)

        if isempty(test_shocks)
            if method == :rotation
                shocks = randn(size(RRR, 2), horizon) # shocks getting drawn are structural shocks
            end
        else
            shocks = test_shocks
        end

        for t = 1:horizon
            out     = vec(X̂' * β) + Σ_chol * shocks[:, t] # X̂ normally would be [X̂ 0 0; 0 X̂ 0; 0 0 X̂] if nobs = 3, but this way of coding it results in less memory storage
            ŷ[t, :] = out

            X̂       = vcat(1., out, X̂[1 + 1:k - nobs]) # XXl = X̂[1 + 1:k - nobs]
        end
    else
        if method == :rotation
            nshocks = size(RRR, 2)
            ŷ       = Array{S, 3}(nshocks, horizon, nobs)
            shocks  = zeros(S, nshocks, horizon)

            for i = 1:nshocks
                shocks[1, i] = flip_shocks ? sqrt(QQ[i, i]) :
                    -sqrt(QQ[i, i]) # a negative 1 s.d. shock by default
                for t = 1:horizon
                    out        = vec(X̂' * β) + Σ_chol * shocks[:, t]
                    ŷ[i, t, :] = out

                    X̂          = vcat(1., out, X̂[1 + 1:k - nobs]) # XXl = X̂[1 + 1:k - nobs]
                end
            end
            ŷ = permutedims(ŷ, [3, 2, 1])
        end

    end

    return ŷ
end
