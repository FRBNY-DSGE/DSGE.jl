"""
```
function impulse_responses(m::AbstractDSGEVARModel{S}, data::AbstractArray{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0 ,use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none)
```
computes the VAR impulse responses identified by the DSGE
```
sₜ = TTT × sₜ₋₁ + RRR × impact[:, i],
yₜ = ZZ × sₜ + DD + MM × impact[:, i],
```
where `impact[:, i]` is a linear combination of
(orthogonal) structural shocks `ϵₜ ∼ 𝒩 (0, I)`, and
`MM × impact[:, i]` are the correlated measurement errors.

We draw a β and Σᵤ from the posterior implied by the DSGE
and data, and we then compute normal VAR impulse responses given those
coefficients and error covariance matrix. The weight placed
on the DSGE is encoded by the field `λ` of the DSGEVAR object `m`.

Given β, Σᵤ, we compute impulse responses to the VAR system
```
ŷₜ₊₁ = X̂ₜ₊₁β + uₜ₊₁,
```
where `X̂ₜ₊₁` are the lags of observables in period `t + 1`, i.e. `yₜ, yₜ₋₁, ..., yₜ₋ₚ`
using one of the available identification methods for VARs

### Inputs
* `method::Symbol`: The available methods are `:cholesky`, `:maxBC`, and `:choleskyLR`.
    See the docstrings `impulse_responses` for VARs specifically.
* `n_obs_shock::Int`: The index of the observable corresponding to the orthogonalized shock
    causing the impulse response.
"""
function impulse_responses(m::AbstractDSGEVARModel{S}, data::AbstractArray{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0, use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}
    β, Σ = compute_system(m, data; verbose = verbose, use_intercept = use_intercept)
    Σ += Σ'
    Σ ./= 2

    return impulse_responses(β, Σ, n_obs_shock,
                             horizon > 0 ? horizon : impulse_response_horizons(m);
                             method = method, use_intercept = use_intercept,
                             flip_shocks = flip_shocks)
end

"""
```
function impulse_responses(m::AbstractDSGEVARModel{S}, data::AbstractArray{S},
    X̂::Matrix{S} = Matrix{S}(undef, 0, 0);
    horizon::Int = 0, MM::Matrix{S} = Matrix{S}(undef, 0, 0),
    flip_shocks::Bool = false, draw_shocks::Bool = false,
    verbose::Symbol = :none) where {S <: Real}
```
computes the VAR impulse responses identified by the DSGE
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

The shock `uₜ₊₁` is identified by assuming
```
Σᵤ = 𝔼[u × u'] = chol(Σᵤ) × Ω × ϵₜ,
```
where the rotation matrix `Ω` is the `Q` matrix from a QR decomposition
of the impact response matrix corresponding to the state space system, i.e.
```
Ω, _ = qr(∂yₜ / ∂ϵₜ').
```

For reference, see Del Negro and Schorfheide (2004),
Del Negro and Schorfheide (2006), and Del Negro and Schorfheide (2009).

### Inputs
* `X̂::Matrix{S}`: covariates for the first "forecast" period
    of the impulse response, i.e. if we have a VAR with `p` lags, then
```
X̂ = [1, ŷₜ, ŷₜ₋₁, ..., ŷₜ₋ₚ₊₁]
```
so that, when β is the vector of VAR coefficients, then
```
𝔼[ŷₜ₊₁] = kron(I, X̂') * β.
```
Internally, we do equivalent matrix operations to avoid allocating
the Kronecker product.

****
NOTE: this function generally involves taking random draws from
probability distributions, so seeds need to be set
to achieve reproducibility.
****

### Keywords
* `horizon::Int`: horizon of impulse responses
* `flip_shocks::Bool`: impulse response shocks are negative by default. Set to `true` for
    a positive signed shock.
* `draw_shocks::Bool`: true if you want to draw shocks along the entire horizon
* `verbose::Symbol`: quantity of output desired
"""
function impulse_responses(m::AbstractDSGEVARModel{S}, data::AbstractArray{S},
                           X̂::Vector{S} = Vector{S}(undef, 0);
                           horizon::Int = 0, MM::Matrix{S} = Matrix{S}(undef, 0, 0),
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
                           verbose::Symbol = :none) where {S <: Real}

    # Prepare X̂
    nobs = size(data, 1)
    k = nobs * get_lags(m) + 1
    if isempty(X̂)
        XX = lag_data(data, get_lags(m); use_intercept = true)
        X̂ = vcat(1, data[:, end], XX[end, 1+1:k - nobs])
    end
    h = (horizon > 0) ? horizon : impulse_response_horizons(m) # grab horizons

    # Compute underlying state space system
    system = compute_system(m; get_system = true, use_intercept = true)

    # Compute VAR coefficients
    β, Σ = compute_system(m, data; verbose = verbose, use_intercept = true)
    Σ += Σ'
    Σ ./= 2.

    # Use rotation identification
    if isempty(MM)
        if hasmethod(measurement_error, (typeof(m),))
            _, MM = measurement_error(m)
        else
            MM = zeros(S, n_observables(m), n_shocks(m))
        end
    end
    return impulse_responses(system[:TTT], system[:RRR], system[:ZZ], system[:DD], MM,
                             system[:QQ], k, β, Σ, X̂, h; flip_shocks = flip_shocks,
                             draw_shocks = draw_shocks)
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
                           X̂::Matrix{S}, horizon::Int;
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
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

For reference, see Del Negro and Schorfheide (2004),
Del Negro and Schorfheide (2006), and Del Negro and Schorfheide (2009).
"""
function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, QQ::Matrix{S},
                           k::Int, β::Matrix{S}, Σ::Matrix{S},
                           X̂::Vector{S}, horizon::Int;
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
                           test_shocks::Matrix{S} =
                           Matrix{S}(undef, 0, 0)) where {S <: Real}

    # Compute impulse responses of predicted values for each β, Σ, and rotation
    nobs = size(ZZ, 1)
    a0_m = convert(Matrix{S},
                   dropdims(impulse_responses(TTT, RRR, ZZ, DD, MM,
                                              sqrt.(QQ), 1, # 1 b/c just want impact and sqrt -> 1 sd shock
                                              accumulate = false); dims = 2)')
    rotation, _ = qr(a0_m)
    Σ_chol = cholesky(Σ).L * rotation' # mapping from structural shocks to innovations in VAR

    if draw_shocks || !isempty(test_shocks)
        ŷ = Matrix{S}(undef, nobs, horizon)

        if isempty(test_shocks)
            shocks = randn(size(RRR, 2), horizon) # shocks getting drawn are structural shocks
        else
            shocks = test_shocks
        end

        for t = 1:horizon
            out     = vec(X̂' * β) + Σ_chol * shocks[:, t] # X̂ normally would be [X̂ 0 0; 0 X̂ 0; 0 0 X̂] if nobs = 3, but this way of coding it results in less memory storage
            ŷ[:, t] = out

            X̂       = vcat(1., out, X̂[1 + 1:k - nobs]) # XXl = X̂[1 + 1:k - nobs]
        end
    else
            nshocks = size(RRR, 2)
            ŷ       = Array{S, 3}(undef, nobs, horizon, nshocks)
            shocks  = zeros(S, nshocks)

            for i = 1:nshocks
                shocks[i] = flip_shocks ? sqrt(QQ[i, i]) :
                    -sqrt(QQ[i, i]) # a negative 1 s.d. shock by default
                out        = vec(X̂' * β) + Σ_chol * shocks # do impact separately
                shocks[i]  = 0. # set back to zero
                ŷ[:, 1, i] = out
                X̂          = vcat(1., out, X̂[1 + 1:k - nobs]) # XXl = X̂[1 + 1:k - nobs]
                for t = 2:horizon
                    out        = vec(X̂' * β)
                    ŷ[:, t, i] = out
                    X̂          = vcat(1., out, X̂[1 + 1:k - nobs]) # XXl = X̂[1 + 1:k - nobs]
                end
            end
    end

    return ŷ
end
