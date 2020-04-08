"""
```
function impulse_responses(m::AbstractDSGEVECMModel{S}, data::AbstractArray{S},
                           coint_mat::AbstractMatrix{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}

function impulse_responses(m::AbstractDSGEVECMModel{S}, coint_mat::AbstractMatrix{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0, use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}
```
computes the VECM impulse responses identified by the DSGE
```
sₜ = TTT × sₜ₋₁ + RRR × impact[:, i],
yₜ = ZZ × sₜ + DD + MM × impact[:, i],
```
where `impact[:, i]` is a linear combination of
(orthogonal) structural shocks `ϵₜ ∼ 𝒩 (0, I)`, and
`MM × impact[:, i]` are the correlated measurement errors.

We draw a β and Σᵤ from the posterior implied by the DSGE
and data, and we then compute normal VECM impulse responses given those
coefficients, innovations variance-covariance matrix, and
the matrix specifying cointegrating relationships in observables. The weight placed
on the DSGE is encoded by the field `λ` of the DSGEVECM object `m`.

Given β, Σᵤ, we compute impulse responses using one of the
available identifiction strategies to the VECM system
```
Δŷₜ₊₁ = eₜβₑ + X̂ₜ₊₁βᵥ + uₜ₊₁,
```
where `βₑ` are the coefficients for the error correction terms;
`eₜ₊₁` are the error correction terms specifying the cointegrating relationships;
`βᵥ` are the coefficients for the VAR terms;
`X̂ₜ₊₁` are the lags of observables in period `t + 1`, i.e. `yₜ, yₜ₋₁, ..., yₜ₋ₚ`,
and `uₜ₊₁ ∼ 𝒩 (0, Σ)`.

If the second function is used (where `data` is not an input), then we assume
the user wants to compute the VECM approximation of the DSGE,
regardless of the `λ` value in `m`. Note that this function will not
update the value of `λ` in `m` (even though we are computing the DSGE-VECM(∞) approximation).

### Inputs
* `coint_mat::AbstractMatrix{S}`: matrix specifying the cointegrating relationships
    in observables. Given a matrix `data` with dimensions `n_observables × T`,
    multiplying `coint_mat * data` should yield a `n_coint × T` matrix, where
    `n_coint` is the number of cointegrating relationships and `T` is
    the number of periods of data.
* `method::Symbol`: The available methods are `:cholesky`, `:maxBC`, and `:choleskyLR`.
    See the docstrings `impulse_responses` for VECMs specifically.
* `n_obs_shock::Int`: The index of the observable corresponding to the orthogonalized shock
    causing the impulse response.

### Keyword Arguments
* `horizon::Int`: the desired horizon of the impulse responses.
* `use_intercept::Bool`: use an intercept term for the VECM approximation
* `flip_shocks::Bool`: default is a "negative" impulse response on impact.
    Set to `true` for the positive impulse response.
"""
function impulse_responses(m::AbstractDSGEVECMModel{S}, data::AbstractArray{S},
                           coint_mat::AbstractMatrix{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = impulse_response_horizons(m),
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}
    β, Σ = compute_system(m, data; verbose = verbose)
    Σ += Σ'
    Σ ./= 2

    return impulse_responses(β, Σ, coint_mat, n_obs_shock, horizon;
                             method = method, flip_shocks = flip_shocks,
                             use_intercept = true)
end


function impulse_responses(m::AbstractDSGEVECMModel{S}, coint_mat::AbstractMatrix{S}, method::Symbol,
                           n_obs_shock::Int; horizon::Int = 0, use_intercept::Bool = false,
                           flip_shocks::Bool = false, verbose::Symbol = :none) where {S <: Real}
    β, Σ = compute_system(m; verbose = verbose, use_intercept = use_intercept)
    Σ += Σ'
    Σ ./= 2

    return impulse_responses(β, Σ, coint_mat, n_obs_shock,
                             horizon > 0 ? horizon : impulse_response_horizons(m);
                             method = method, use_intercept = use_intercept,
                             flip_shocks = flip_shocks)
end

"""
```
function impulse_responses(m::AbstractDSGEVECMModel{S}, data::AbstractArray{S},
    X̂::Matrix{S} = Matrix{S}(undef, 0, 0);
    horizon::Int = 0, MM::Matrix{S} = Matrix{S}(undef, 0, 0),
    flip_shocks::Bool = false, draw_shocks::Bool = false,
    verbose::Symbol = :none) where {S <: Real}
```
computes the VECM impulse responses identified by the DSGE
```
sₜ = TTT × sₜ₋₁ + RRR × impact[:, i],
yₜ = ZZ × sₜ + DD + MM × impact[:, i],
```
where `impact[:, i]` is a linear combination of
(orthogonal) structural shocks `ϵₜ ∼ 𝒩 (0, I)`, and
`MM × impact[:, i]` are the correlated measurement errors.

The VECM impulse responses are computed according to
```
Δŷₜ₊₁ = eₜβₑ + X̂ₜ₊₁βᵥ + uₜ₊₁,
```
where `βₑ` are the coefficients for the error correction terms;
`eₜ` are the error correction terms specifying the cointegrating relationships;
`βᵥ` are the coefficients for the VAR terms;
`X̂ₜ₊₁` are the lags of observables in period `t + 1`, i.e. `yₜ, yₜ₋₁, ..., yₜ₋ₚ`,
and `uₜ₊₁ ∼ 𝒩 (0, Σ)`. Note these impulses responses are *not*
computed in deviations from the baseline forecast `Δŷₜ₊₁ = eₜ₊₁βₑ + X̂ₜ₊₁βᵥ`.
To compute these impulse responses, use the keyword `deviations`.

The shock `uₜ₊₁` is identified by assuming
```
Σᵤ = 𝔼[u × u'] = chol(Σᵤ) × Ω × ϵₜ,
```
where the rotation matrix `Ω` is the `Q` matrix from a QR decomposition
of the impact response matrix corresponding to the state space system, i.e.
```
Ω, _ = qr(∂yₜ / ∂ϵₜ').
```
The impact response matrix is constructed using only the stationary component of the
state space system and ignores the cointegration components of `ZZ` and `DD`.

For reference, see Del Negro and Schorfheide (2004),
Del Negro and Schorfheide (2006), and Del Negro and Schorfheide (2009).

### Inputs
* `coint_mat::Matrix{S}`: matrix specifying the cointegrating relationships
    in the actual `data` matrix. Evaluating `coint_mat * data` should yield
    a time series of the cointegrating relationships.
* `X̂::Matrix{S}`: covariates for the first "forecast" period
    of the impulse response, i.e. if we have a VECM with `p` lags, then
```
X̂ = [eₜ, 1, ŷₜ, ŷₜ₋₁, ..., ŷₜ₋ₚ₊₁]
```
so that, when β is the vector of VECM coefficients, then
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
* `deviations::Bool`: set true to compute the impulse response in deviations
    rather than as a forecast. Mechnically, we ignore `X̂` (treated as zeros)
    and the intercept term.
* `verbose::Symbol`: quantity of output desired
"""
function impulse_responses(m::AbstractDSGEVECMModel{S}, data::AbstractArray{S},
                           coint_mat::Matrix{S}, X̂::Vector{S} = Vector{S}(undef, 0);
                           horizon::Int = 0, MM::Matrix{S} = Matrix{S}(undef, 0, 0),
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
                           deviations::Bool = false, verbose::Symbol = :none) where {S <: Real}

    # Prepare X̂
    n_obs = size(data, 1)
    n_coint = size(coint_mat, 1)
    k = n_coint + n_obs * get_lags(m) + 1
    if isempty(X̂)
        XX = lag_data(data, get_lags(m); use_intercept = true)
        addcoint = coint_mat * data[:, end]
        X̂ = vcat(addcoint, 1, data[:, end], XX[end, 1+1:k - n_obs])
    end
    h = (horizon > 0) ? horizon : impulse_response_horizons(m) # grab horizons

    # Compute underlying state space system
    system = compute_system(m; get_system = true, use_intercept = true)

    # Compute VECM coefficients
    β, Σ = compute_system(m, data; verbose = verbose)
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
                             system[:QQ], k, n_obs, n_coint, β, Σ,
                             coint_mat, X̂, h; flip_shocks = flip_shocks,
                             draw_shocks = draw_shocks, deviations = deviations)
end

"""
```
function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, QQ::Matrix{S},
                           k::Int, n_obs::Int, n_coint::Int, β::Matrix{S}, Σ::Matrix{S},
                           coint_mat::Matrix{S}, horizon::Int, X̂::Matrix{S} = zeros(S, k);
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
                           deviations::Bool = false,
                           test_shocks::Matrix{S} =
                           Matrix{S}(undef, 0, 0)) where {S<:Real}
```
computes the VECM impulse responses identified by the state space system
```
sₜ = TTT × sₜ₋₁ + RRR × ϵₜ
yₜ = ZZ × sₜ + DD + MM × ϵₜ
```
where `ϵₜ ∼ 𝒩 (0, QQ)` and `MM × ϵₜ` are the correlated measurement errors.

Consider the VECM
```
Δŷₜ₊₁ = eₜβₑ + X̂ₜ₊₁βᵥ + uₜ₊₁,
```
where `βₑ` are the coefficients for the error correction terms;
`eₜ` are the error correction terms specifying the cointegrating relationships;
`βᵥ` are the coefficients for the VAR terms (including the intecept)o;
`X̂ₜ₊₁` are the lags of observables in period `t + 1`, i.e. `yₜ, yₜ₋₁, ..., yₜ₋ₚ₊₁`;
and `uₜ₊₁ ∼ 𝒩 (0, Σ)`. Note these impulses responses are *not*
computed in deviations from the baseline forecast `Δŷₜ₊₁ = eₜ₊₁βₑ + X̂ₜ₊₁βᵥ`.
To compute these impulse responses, set the keyword `deviations = true`.

The shock `uₜ₊₁` is identified via
```
Σᵤ = 𝔼[u × u'] = chol(Σᵤ) × Ω × ϵₜ,
```
where the rotation matrix `Ω` is the `Q` matrix from a QR decomposition
of the impact response matrix corresponding to the state space system, i.e.
```
Ω, _ = qr(∂yₜ / ∂ϵₜ').
```
The impact response matrix is constructed using only the stationary component of the
state space system and ignores the cointegration components of `ZZ` and `DD`.

The data are assumed to have dimensions `n_obs × T`, and
the cointegration relationships in the data are given by `coint_mat * data`, where
`coint_mat` has dimensions `n_coint × n_obs`. The variable `k` is the
number of total regressors in the VECM, including cointegration terms.

For reference, see Del Negro and Schorfheide (2004), Del Negro and Schorfheide (2006),
and Del Negro, Schorfheide, Smets, and Wouters (2007).
"""
function impulse_responses(TTT::Matrix{S}, RRR::Matrix{S}, ZZ::Matrix{S},
                           DD::Vector{S}, MM::Matrix{S}, QQ::Matrix{S},
                           k::Int, n_obs::Int, n_coint::Int, β::Matrix{S}, Σ::Matrix{S},
                           coint_mat::Matrix{S}, horizon::Int, X̂::Vector{S} = zeros(S, k);
                           flip_shocks::Bool = false, draw_shocks::Bool = false,
                           deviations::Bool = false,
                           test_shocks::Matrix{S} =
                           Matrix{S}(undef, 0, 0)) where {S <: Real}

    # Grab stationary components
    ZZ = ZZ[1:n_obs, :]
    DD = DD[1:n_obs]
    MM = MM[1:n_obs, :]

    # Compute impulse responses of predicted values for each β, Σ, and rotation
    a0_m = convert(Matrix{S},
                   dropdims(impulse_responses(TTT, RRR, ZZ, DD, MM,
                                              sqrt.(QQ), 1, # 1 b/c just want impact and sqrt -> 1 sd shock
                                              accumulate = false); dims = 2)')
    rotation, r_a = qr(a0_m)
    rotation = sign.(diag(r_a)) .* rotation' # normalizes each row (change signs) so that lower diagonal (r_a') has positive diagonal elements
    Σ_chol = cholesky(Σ).L * rotation # mapping from structural shocks to innovations in VECM

    # Remove intercept term if in deviations
    if deviations
        β = β[vcat(1:n_coint, n_coint + 2:size(β, 1)), :]
        X̂ = zeros(S, size(β, 1))
        k -= 1
    end

    if draw_shocks || !isempty(test_shocks)
        ŷ = Matrix{S}(undef, n_obs, horizon)

        if isempty(test_shocks)
            shocks = randn(size(RRR, 2), horizon) # shocks getting drawn are structural shocks
            if flip_shocks
                @warn "Shocks are being drawn, so flipping shocks does not make sense in this context. Ignoring `flip_shocks` . . ."
            end
        else
            shocks = test_shocks
        end
        for t = 1:horizon
            out      = vec(X̂' * β) + Σ_chol * shocks[:, t] # X̂ normally would be [X̂ 0 0; 0 X̂ 0; 0 0 X̂] if n_obs = 3,
            ŷ[:, t]  = out                                 # but this way of coding it results in less memory storage
            addcoint = X̂[1:n_coint] + coint_mat * out      # Predicted cointegration terms

            X̂ = deviations ? vcat(addcoint,  out, X̂[n_coint + 1:k - n_obs]) :
                vcat(addcoint, 1.,  out, X̂[n_coint + 1 + 1:k - n_obs]) # XXl = X̂[n_coint + 1 + 1:k - n_obs]
        end
    else
        nshocks = size(RRR, 2)
        ŷ       = Array{S, 3}(undef, n_obs, horizon, nshocks)
        old_X̂   = X̂
        shocks  = zeros(S, nshocks)

        for i = 1:nshocks
            X̂ = copy(old_X̂)
            shocks[i] = flip_shocks ? sqrt(QQ[i, i]) :
                -sqrt(QQ[i, i]) # a negative 1 s.d. shock by default
            out        = vec(X̂' * β) + Σ_chol * shocks # do impact separately
            shocks[i]  = 0. # set back to zero
            ŷ[:, 1, i] = out
            addcoint   = X̂[1:n_coint] + coint_mat * out
            X̂          = deviations ? vcat(addcoint,  out, X̂[n_coint + 1:k - n_obs]) :
                vcat(addcoint, 1.,  out, X̂[n_coint + 1 + 1:k - n_obs]) # XXl = X̂[n_coint + 1 + 1:k - n_obs]
            for t = 2:horizon
                out        = vec(X̂' * β)
                ŷ[:, t, i] = out
                addcoint   = X̂[1:n_coint] + coint_mat * out
                X̂          = deviations ? vcat(addcoint,  out, X̂[n_coint + 1:k - n_obs]) :
                    vcat(addcoint, 1.,  out, X̂[n_coint + 1 + 1:k - n_obs]) # XXl = X̂[n_coint + 1 + 1:k - n_obs]
            end
        end
    end

    return ŷ
end
