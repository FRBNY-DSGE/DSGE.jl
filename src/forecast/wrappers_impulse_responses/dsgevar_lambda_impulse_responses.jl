"""
```
function impulse_responses(m::AbstractDSGEVARModel{S}, paras::Matrix{S},
                           data::Matrix{S}, input_type::Symbol, method::Symbol;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           n_obs_shock::Int = 1, draw_shocks::Bool = false,
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false, test_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
```
computes the VAR impulse responses identified by the state space system
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
coefficients and innovations variance-covariance matrix.

****
NOTE: this function generally involves taking random draws from
probability distributions, so seeds need to be set
to achieve reproducibility.
****

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
* `n_obs_shock::Int`: the index of the observable corresponding to the orthogonalized shock causing the impulse response.

### Keywords
* `parallel::Bool`: use parallel workers or not
* `frequency_band::Tuple{S,S}`: See `?maxBC_shock`.
* `n_obs_shock::Int`: Index of observable to be shocked when using a Cholesky-based impulse response
* `draw_shocks::Bool`: true if you want to draw shocks along the entire horizon
* `flip_shocks::Bool`: impulse response shocks are negative by default. Set to `true` for
    a positive signed shock.
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
                           n_obs_shock::Int = 1, draw_shocks::Bool = false,
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false, test_meansbands::Bool = false,
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
    nobs = size(data, 1)
    k = nobs * lags + 1
    h = impulse_response_horizons(m)

    # Preompute X̂ and MM for the rotation identification
    if method == :rotation
        XX = lag_data(data, lags; use_intercept = true)
        X̂ = vcat(1, data[:, end], XX[end, 1+1:k - nobs])

        # Get measurement error
        if hasmethod(measurement_error, (typeof(m),))
            _, MM = measurement_error(m)
        else
            MM = zeros(S, n_observables(m), n_shocks(m))
        end
    end

    dsgevarrotationirf_method = if method == :rotation
        function _dsgevar_λ_rotation_irf_(para)
            DSGE.update!(m, para)
            return impulse_responses(m, data, X̂; horizon = h, MM = MM,
                                     flip_shocks = flip_shocks, draw_shocks = draw_shocks,
                                     verbose = verbose)
        end
    elseif method in [:cholesky, :cholesky_long_run, :choleskyLR, :maxBC,
                      :maximum_business_cycle_variance]
        function _dsgevar_λ_irf_(para)
            DSGE.update!(m, para)
            return impulse_responses(m, data, method, n_obs_shock; horizon = h,
                                     use_intercept = true, flip_shocks = flip_shocks)
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
        if method == :rotation && !draw_shocks
            for (shock, shock_i) in get_shocks(m)
                for (name, name_i) in get_observables(m)
                    # irf_output is Vector{nperiod x nobs x nshocks} -> for each observable,
                    # we want to select its specific IRF, i.e. map(x -> x[:,obs_index, shock_index]).
                    # This creates a nperiod x ndraws matrix, which we want to transpose
                    # to get a ndraws x nperiod matrix
                    single_var = Matrix(reduce(hcat, map(x -> x[:, name_i, shock_i], irf_output))')
                    means[!, Symbol(name, :__, shock)] = vec(mean(single_var, dims = 1))
                    bands[Symbol(name, :__, shock)]    = find_density_bands(single_var, density_bands;
                                                                            minimize = minimize)
                end
            end
        else
            for (name, name_i) in get_observables(m)
                # irf_output is Vector{nperiod x nobs} -> for each observable,
                # we want to select its specific IRF, i.e. map(x -> x[:,obs_index]).
                # This creates a nperiod x ndraws matrix, which we want to transpose
                # to get a ndraws x nperiod matrix
                single_var = Matrix(reduce(hcat, map(x -> x[:, name_i], irf_output))')
                means[!, name] = vec(mean(single_var, dims = 1))
                bands[name]    = find_density_bands(single_var, density_bands;
                                                    minimize = minimize)
            end
        end
        mb = MeansBands(metadata, means, bands)

        # Save MeansBands
        if !test_meansbands
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
        end
        return mb
    else
        if method in [:rotation]
            # Reshape irf_output to nobs x nperiod x nshock x ndraw
            return cat(irf_output..., dims = ndims(irf_output[1]) + 1)
        else
            # Reshape irf_output to nobs x nperiod x ndraw
            return cat(irf_output..., dims = 3)
        end
    end
end

function impulse_responses(m::AbstractDSGEVARModel{S}, paras::Vector{S},
                           data::Matrix{S}, input_type::Symbol, method::Symbol;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           n_obs_shock::Int = 1, draw_shocks::Bool = false,
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    return impulse_responses(m, reshape(paras, 1, length(paras)),
                             data, input_type, method; parallel = parallel,
                             frequency_band = frequency_band, n_obs_shock = n_obs_shock,
                             draw_shocks = draw_shocks, flip_shocks = flip_shocks,
                             density_bands = density_bands,
                             create_meansbands = create_meansbands, minimize = minimize,
                             forecast_string = forecast_string, verbose = verbose)
end
