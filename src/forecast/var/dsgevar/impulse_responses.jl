"""
```
function impulse_responses(m, input_type, method,
                           lags, observables, shocks,
                           n_obs_var; parallel = false,
                           frequency_band = (2*π/32, 2*π/6),
                           flip_shocks = false,
                           density_bands = [.5, .6, .7, .8, .9],
                           compute_meansbands = false,
                           minimize = true,
                           forecast_string = "",
                           verbose = :high) where {S<:Real}
function impulse_responses(m, paras, input_type, method,
                           lags, observables, shocks,
                           n_obs_var; parallel = false,
                           frequency_band = (2*π/32, 2*π/6),
                           flip_shocks = false,
                           density_bands = [.5, .6, .7, .8, .9],
                           compute_meansbands = false,
                           minimize = true,
                           forecast_string = "",
                           verbose = :high) where {S<:Real}
```
computes the impulse responses of a VAR(p) approximation to a DSGE.

### Inputs
* `m::AbstractDSGEModel`: DSGE model object
* `paras::Matrix{S}` or `paras::Vector{S}`: parameters to calibrate the DSGE
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
    These shocks must be in `m`
* `n_obs_var::Int`: the index of the observable to be shocked by
    the reduced-form impulse response to the VAR system.

### Keywords
* `parallel::Bool`: use parallel workers or not
* `frequency_band::Tuple{S,S}`: See `?maxBC_shock`.
* `flip_shocks::Bool`: impulse response shocks are negative by default. Set to `true` for
    a positive signed shock.
* `density_bands::Vector{Float64}`: bands for full-distribution IRF computations
* `compute_meansbands::Bool`: set to `true` to save output as a `MeansBands` object.
* `minimize::Bool`: choose shortest interval if true, otherwise just chop off lowest and
    highst (percent/2)
* `forecast_string::String`: string tag for identifying this impulse response
* `verbose::Symbol`: quantity of output desired

"""
function impulse_responses(m::AbstractDSGEModel, input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_var::Int;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           compute_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    return impulse_responses(m, load_draws(m, input_type),
                             input_type, method, lags, observables,
                             shocks, n_obs_var;
                             parallel = parallel,
                             frequency_band = frequency_band,
                             flip_shocks = flip_shocks,
                             density_bands = density_bands,
                             compute_meansbands = compute_meansbands,
                             minimize = minimize,
                             forecast_string = forecast_string,
                             verbose = verbose) where {S<:Real}
end

function impulse_responses(m::AbstractDSGEModel, paras::Vector{S},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_var::Int;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           compute_meansbands::Bool = false,
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
                             compute_meansbands = compute_meansbands,
                             minimize = minimize,
                             forecast_string = forecast_string,
                             verbose = verbose) where {S<:Real}
end

function impulse_responses(m::AbstractDSGEModel, paras::Matrix{S},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           shocks::Vector{Symbol},
                           n_obs_var::Int;
                           parallel::Bool = false,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           compute_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           verbose::Symbol = :high) where {S<:Real}
    if n_obs_var <= 0
        error("To use method $method, user must specify the index of" *
              " the target observable with keyword n_obs_var.")
    end

    # Set up computation method
    mapfcn = parallel ? pmap : map
    h = impulse_response_horizons(m)

    # Compute VAR coefficients implied by DSGE
    paras = mapslices(x -> [vec(x)], paras, dims = 2)
    var_output =
        mapfcn(para -> dsge_var!(m, para, observables, shocks, lags), paras)

    # Reformat output
    β_draws = map(x -> x[1], var_output)
    Σ_draws = map(x -> x[2], var_output)

    # Compute IRFs
    irf_output =
        mapfcn((β, Σ) ->
               impulse_responses(β, Σ, n_obs_var, h; method = method,
                                 include_constant = false,
                                 flip_shocks = flip_shocks),
               β_draws, Σ_draws)

    if compute_meansbands
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
            :_cholesky
        elseif method == :maxBC || method == :maximum_business_cycle_variance
            :_maxBC
        else
            :_choleskyLR
        end

        var_names = Symbol(join(string.(observables), "_"))
        fp = get_meansbands_output_file(m, input_type, :none,
                                        Symbol(:dsgevarirf, :obs_,
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


# """
# ```
# impulse_responses(β, Σ, n_obs_shock, horizon, shock_size = 1;
#     method = :cholesky, flip_shocks = false, include_constant = true,
#     frequency_band = (2π/32, 2π/6)) where {S<:Real}
# ```
# computes the impulse responses of a VAR system represented in the form

# ```
# yₜ = Xₜβ + ϵₜ,
# ```
# where `Xₜ` stacks the lags of yₜ (with dimensions n_observables x n_regressors), and

# ```
# ϵₜ ∼ 𝒩 (0, Σ).
# ```

# ### Inputs
# * `β::AbstractMatrix{S}`: coefficient matrix
# * `Σ::AbstractMatrix{S}`: innovations covariance matrix
# * `n_obs_shock::Int`: index of the observable to be shocked
# * `shock_size::S`: number of standard deviations of the shock

# ### Keywords
# * `method::Symbol`: type of impulse response to compute. The available options are
#     `:cholesky` (default), `:maximum_business_cycle_variance` or `:maxBC`, and
#     `:cholesky_long_run` or `:choleskyLR`. See `?cholesky_shock`, `?maxBC_shock`,
#     and `?cholesky_long_run_shock`.
# * `flip_shocks::Bool`: by default, we compute the impulse responses to a negative shock.
#     Set `flip_shocks = true` to obtain a positive shock.
# * `include_constant::Bool`: `impulse_responses` assumes `β` has constant term(s). If there
#     are no such terms, then `include_constant` must be set to `false`.
# * `frequency_band::Tuple{S,S}`: See `?maxBC_shock`.

# ### Outputs
# * `Y::AbstractMatrix`: Impulse response matrix with dimensions horizons x n_observables
# """
# function impulse_responses(β::AbstractMatrix{S}, Σ::AbstractMatrix{S}, n_obs_shock::Int,
#                            horizon::Int, shock_size::S = one(S);
#                            method::Symbol = :cholesky,
#                            flip_shocks::Bool = false,
#                            include_constant::Bool = true,
#                            frequency_band::Tuple{S,S} =
#                            (2*π/32, 2*π/6)) where {S<:Real}

#     # Compute dimensions
#     n = size(β,2)
#     lags = convert(Int, include_constant ? (size(β,1) - 1) / n : size(β,1) / n)

#     # Compute impact based on IRF type
#     Y = zeros(lags + horizon, n)
#     Y[lags + 1, :] = if method == :cholesky
#         cholesky_shock(Σ, n, n_obs_shock, shock_size;
#                        flip_shocks = flip_shocks)
#     elseif method == :maximum_business_cycle_variance || method == :maxBC
#         maxBC_shock(β, Σ, n, n_obs_shock, shock_size, lags, frequency_band;
#                     flip_shocks = flip_shocks)
#     elseif method == :choleskyLR || method == :cholesky_long_run
#         cholesky_long_run_shock(β, Σ, n_obs_shock, n, lags, shock_size;
#                                 flip_shocks = flip_shocks)
#     else
#         error("IRF method $(string(method)) has not been implemented.")
#     end

#     # For efficiency
#     if include_constant
#         β = @views β[2:end, :]
#     end

#     # Compute impulse response
#     for t = 2:horizon
#         xT = reshape(Y[lags + t - 1:-1:lags + t - lags, :]', lags * n, 1)'
#         Y[lags + t, :] = xT * β
#     end

#     return Y[lags + 1:end, :]
# end

# """
# ```
# cholesky_shock(Σ, n, n_obs_shock, shock_size, flip_shocks = false) where {S<:Real}
# ```
# computes a Cholesky-identified shock to the specified observable.

# ### Inputs
# * `Σ::AbstractMatrix{S}`: innovations covariance matrix
# * `n::Int`: number of observables
# * `n_obs_shock::Int`: index of the observable to be shocked
# * `shock_size::S`: number of standard deviations of the shock

# ### Keywords
# * `flip_shocks::Bool`: by default, we compute the impulse responses to a negative shock.
#     Set `flip_shocks = true` to obtain a positive shock.
# """
# function cholesky_shock(Σ::Matrix{S}, n::Int, n_obs_shock::Int,
#                         shock_size::S; flip_shocks::Bool = false) where {S<:Real}
#     cholmat = cholesky((Σ + Σ') ./ 2).L
#     vec_shock = zeros(n)
#     vec_shock[n_obs_shock] = flip_shocks ? shock_size : -shock_size # negative by DSGE convention
#     return (cholmat * vec_shock)'
# end

# """
# ```
# maxBC_shock(β, Σ, n, n_obs_shock, shock_size, lags, frequency_band,
#     flip_shocks = false) where {S<:Real}
# ```
# maximizes the business cycle variance explained by the observable
# whose index is specified by `n_obs_shock` and between the
# frequencies specified by `frequency_band`.

# ### Inputs
# * `β::AbstractMatrix{S}`: coefficient matrix
# * `Σ::AbstractMatrix{S}`: innovations covariance matrix
# * `n::Int`: number of observables
# * `n_obs_shock::Int`: index of the observable to be shocked
# * `shock_size::S`: number of standard deviations of the shock
# * `lags::Int`: number of lags in VAR system
# * `frequency_band::Tuple{S,S}`: the frequencies between which the variance of
#     the observable specified by `n_obs_shock` will be maximized.

# ### Keywords
# * `flip_shocks::Bool`: by default, we compute the impulse responses to a negative shock.
#     Set `flip_shocks = true` to obtain a positive shock.
# """
# function maxBC_shock(β::Matrix{S}, Σ::Matrix{S}, n::Int, n_obs_shock::Int, shock_size::S,
#                      lags::Int, frequency_band::Tuple{S,S};
#                      flip_shocks::Bool = false) where {S<:Real}
#     if lags * n < size(β,1)
#         β = @views β[2:end, :]
#     end

#     cholmat = cholesky((Σ + Σ') ./ 2).L
#     increment = abs(frequency_band[1] - frequency_band[2]) / 200.
#     V = zeros(S, n, n) # variance
#     eminusif = zeros(Complex{S}, 1, 1, lags)
#     for f = frequency_band[1]:increment:round(frequency_band[2], digits=10) # not rounding sometimes leads to one fewer loop than desired
#         eminusif[1, 1, :] = exp.(-im .* f .* collect(1:lags))
#         sumB = dropdims(sum(reshape(β', n, n, lags) .*
#                            repeat(eminusif, n, n, 1); dims = 3), dims = 3)
#         invA = (Matrix{Complex{S}}(I, n, n) - sumB) \ cholmat
#         V += reshape(real.(kron(conj(invA[n_obs_shock, :]), invA[n_obs_shock, :])), n, n) .*
#             increment ./ abs(frequency_band[1] - frequency_band[2])
#     end
#     eigout = eigen(V)
#     q = eigout.vectors[:, argmax(eigout.values)]
#     q .*= sign(q[n_obs_shock])
#     q .*= flip_shocks ? shock_size : -shock_size # negative by DSGE convention

#     return (cholmat * q)'
# end

# """
# ```
# cholesky_long_run_shock(β, Σ, n, n_obs_shock, shock_size, lags, frequency_band,
#     flip_shocks = false) where {S<:Real}
# ```
# computes the long-run Cholesky-identified shock to the observable
# specified by `n_obs_shock`.

# Given a VAR system
# ```
# yₜ = B₁yₜ₋₁ + ... + B₁yₜ₋ₚ + Γϵₜ,      ϵₜ ∼ 𝒩 (0, Σ),
# ```

# the long-run covariance matrix is
# ```
# S̃ = B̃⁻¹ Σ (B̃⁻¹)'
# ```

# and the Cholesky identification is given by
# ```
# ΓΓ' = Σ ⇒ Γ = B̃ * cholesky(S̃).
# ```

# ### Inputs
# * `β::AbstractMatrix{S}`: coefficient matrix
# * `Σ::AbstractMatrix{S}`: innovations covariance matrix
# * `n::Int`: number of observables
# * `n_obs_shock::Int`: index of the observable to be shocked
# * `shock_size::S`: number of standard deviations of the shock
# * `lags::Int`: number of lags in VAR system
# * `frequency_band::Tuple{S,S}`: the frequencies between which the variance of
#     the observable specified by `n_obs_shock` will be maximized.

# ### Keywords
# * `flip_shocks::Bool`: by default, we compute the impulse responses to a negative shock.
#     Set `flip_shocks = true` to obtain a positive shock.
# """
# function cholesky_long_run_shock(β::Matrix{S}, Σ::Matrix{S}, n_obs_shock::Int, n::Int,
#                                  lags::Int, shock_size::S;
#                                  flip_shocks::Bool = false) where {S<:Real}
#     if n * lags < size(β, 1)
#         β = β[2:end,:] # don't need the constant
#     end

#     # Compute decomposition
#     B̃ = Matrix{S}(I, n, n) - dropdims(sum(reshape(β', n, n, lags), dims = 3), dims = 3)
#     S̃ = B̃ \ (Σ * inv(B̃)')             # LR covariance = S̃ = B̃⁻¹ * Σ * B̃⁻¹' =>
#     Γ = B̃ * cholesky((S̃ + S̃') ./ 2).L # S = B̃ \ (Σ * B̃⁻¹')

#     # Compute shock
#     vec_shock = zeros(n)
#     vec_shock[n_obs_shock] = flip_shocks ? shock_size : -shock_size # negative by DSGE convention
#     return (Γ * vec_shock)'
# end

function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, impact::Matrix{S}, horizon::Int;
                           accumulate::Bool = false,
                           cum_inds::Union{Int,UnitRange{Int},Vector{Int}} = 0) where {S<:Real}
    # Get dimensions
    nobs   = size(impact, 1)
    nstate = size(TTT, 1)

    # Compute impulse response to identified impact matrix
    irf_results = Matrix{S}(undef, horizon, nobs^2)
    for i = 1:nobs
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
            obs[cum_inds, :] = cumsum(obs[cum_inds, :], dims = 2)
        end
        irf_results[:, 1 + (i - 1) * nobs:i * nobs] = obs'
    end

    return irf_results
end

function impulse_responses_rotation(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                                    DD::Vector{S}, MM::Matrix{S}, impact::Matrix{S},
                                    k::Int, β::Matrix{S}, Σ::Matrix{S},
                                    x̂::Matrix{S}, horizon::Int;
                                    accumulate::Bool = false,
                                    cum_inds::Union{Int,UnitRange{Int},Vector{Int}}
                                    = 0, test_shocks::Matrix{S} =
                                    Matrix{S}(undef, 0, 0)) where {S<:Real}
    aairf = impulse_responses(TTT, RRR, ZZ, DD, MM, impact, 1, # 1 b/c just want impact
                              accumulate = accumulate, cum_inds = cum_inds)
    nobs = size(ZZ, 1)
    a0_m = Matrix{S}(undef, nobs, nobs)
    for i = 1:nobs
        a0_m[i, :] = aairf[1, nobs * (i - 1) + 1:nobs * i]
    end
    β_rotation, _ = qr(a0_m)
    β_rotation    = convert(Matrix{S}, β_rotation')

    # Compute distribution of predicted values for each β, Σ, and rotation
    Σ_chol = cholesky(Σ).L * β_rotation
    ŷ      = Matrix{S}(undef, horizon, nobs)
    shocks = isempty(test_shocks) ? randn(horizon, nobs) : test_shocks

    for t = 1:horizon
        out     = (β * x̂' + Σ_chol * shocks[t,:])'
        ŷ[t, :] = out
        xxl = x̂[1 + 1:k - nobs]
        x̂ = hcat(1., out, reshape(xxl, 1, length(xxl)))
    end

    return ŷ
end
