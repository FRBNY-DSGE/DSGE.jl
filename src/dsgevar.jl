"""
```
dsge_var(m, observables, shocks, lags; zero_DD = false,
    EE = [], MM = zeros(length(observables), length(shocks)),
    get_VAR = true) where {S<:Real}
dsge_var!(m, para, observables, shocks, lags; zero_DD = false,
    EE = [], MM = zeros(length(observables), length(shocks)),
    get_VAR = true) where {S<:Real}
```
are wrappers for computing the VAR(p) approximation of DSGE model `m`.

See `?var_approx_state_space` for the routine which performs the computation.

### Inputs
* `m::AbstractDSGEModel`: model object
* `observables::Vector{Symbol}`: requested observables in VAR. These can be
    `Observables` or `PseudoObservables`, as long as they are in
    the dictionary of observables or pseudo-observables of model object `m`.
* `shocks::Vector{Symbol}`: requested (structural) exogenous shocks. These must be
    located in the dictionary of exogenous shocks of the model object `m`.
* `lags::Int`: number of lags for the VAR approximation.

### Keyword Arguments
* `EE::AbstractMatrix{S}`: measurement error covariance matrix EE
* `MM::AbstractMatrix{S}`: implements correlation between measurement error
    and the exogenous shocks. See `?var_approx_state_space`

### Outputs
If `get_VAR = true`:
* `β`: VAR(p) coefficients
* `Σ`: innovations covariance matrix for the VAR(p) representation
```
yₜ = Xₜβ + μₜ
```
where `Xₜ` appropriately stacks the `p` lags of `yₜ` and `μₜ ∼ 𝒩 (0, Σ)`.

If `get_VAR = false`: we return the limit cross product matrices.
* `yyyyd`: 𝔼[y,y]
* `XXXXd`: 𝔼[y,X(lag rr)]
* `XXyyd`: 𝔼[X(lag rr),X(lag ll)]

Using these matrices, the VAR(p) representation is given by
```
β = XXXXd \ XXyyd
Σ = yyyyd - XXyyd' * β
```
"""
function dsge_var(m::AbstractDSGEModel, observables::Vector{Symbol},
                  shocks::Vector{Symbol}, lags::Int; zero_DD::Bool = false,
                  MM::AbstractMatrix{S} =
                  zeros(S, length(observables), length(shocks)),
                  EE::AbstractMatrix{S} =
                  zeros(S, length(observables), length(observables)),
                  get_VAR::Bool = true) where {S<:Real}
    para = map(x -> x.value, m.parameters)
    return dsge_var!(m, para, observables, shocks, lags;
                     zero_DD = zero_DD, MM = MM, EE = EE, get_VAR = get_VAR)
end

function dsge_var!(m::AbstractDSGEModel, para::Vector{S},
                   observables::Vector{Symbol},
                   shocks::Vector{Symbol}, lags::Int;
                   zero_DD::Bool = false,
                   MM::AbstractMatrix{S} =
                   zeros(S, length(observables), length(shocks)),
                   EE::AbstractMatrix{S} =
                   Matrix{S}(undef,0,0),
                   get_VAR::Bool = true) {S<:Real}
    DSGE.update!(m, para)
    system = compute_system(m)
    system = compute_system(m, system; observables = observables,
                            shocks = shocks, zero_DD = zero_DD)
    if !isempty(EE)
        system.measurement.EE = EE
    end

    return var_approx_state_space(system[:TTT], system[:RRR], system[:QQ], system[:DD],
                                  system[:ZZ], system[:EE], MM, lags; get_VAR = get_VAR)
end

"""
```
var_approx_state_space(TTT, RRR, QQQ, DD, ZZ, EE, MM, p; get_VAR = true) where {S<:Real}
```
computes the VAR(p) approximation of the linear state space system

```
sₜ = TTT * sₜ₋₁ + RRR * ϵₜ,
yₜ = ZZ * sₜ + DD + uₜ,
```
where the disturbances are assumed to follow
```
ϵₜ ∼ 𝒩 (0, QQ),
uₜ = ηₜ + MM * ηₜ,
ηₜ ∼ 𝒩 (0, EE).
```
The `MM` matrix implies
```
cov(ϵₜ, uₜ) = QQ * MM'.
```

### Outputs
If `get_VAR = true`:
* `β`: VAR(p) coefficients
* `Σ`: innovations covariance matrix for the VAR(p) representation
```
yₜ = Xₜβ + μₜ
```
where `Xₜ` appropriately stacks the `p` lags of `yₜ` and `μₜ ∼ 𝒩 (0, Σ)`.

If `get_VAR = false`: we return the limit cross product matrices.
* `yyyyd`: 𝔼[y,y]
* `XXXXd`: 𝔼[y,X(lag rr)]
* `XXyyd`: 𝔼[X(lag rr),X(lag ll)]

Using these matrices, the VAR(p) representation is given by
```
β = XXXXd \ XXyyd
Σ = yyyyd - XXyyd' * β
```
"""
function var_approx_state_space(TTT::AbstractMatrix{S}, RRR::AbstractMatrix{S},
                                QQ::AbstractMatrix{S}, DD::AbstractVector{S},
                                ZZ::AbstractMatrix{S}, EE::AbstractMatrix{S},
                                MM::AbstractMatrix{S},
                                p::Int; get_VAR::Bool = true) where {S<:Real}

    nobs = size(ZZ,1)

    yyyyd = zeros(S, nobs, nobs)
    XXyyd = zeros(S, p * nobs, nobs)
    XXXXd = zeros(S, p * nobs, p * nobs)

    HH = EE + MM * QQ * MM';
    VV = QQ * MM';

    ## Compute p autocovariances

    ## Initialize Autocovariances
    GAMM0 = zeros(S, nobs ^ 2, p + 1)

    GA0 =  solve_discrete_lyapunov(TTT, RRR * QQ * RRR')
    Gl   = ZZ * GA0 * ZZ' + ZZ * RRR * VV + (ZZ * RRR * VV)' + HH
    GAMM0[:, 1] = vec(Gl)

    TTl = copy(TTT)
    GA0ZZ = GA0 * ZZ'
    RRRVV = RRR * VV
    for l = 1:p
        Gl = ZZ * TTl * GA0ZZ + ZZ * TTl * RRRVV # ZZ * (TTl * GA0Z) * ZZ' + ZZ * (TTl * RRR * VV)
        GAMM0[:, l+1] = vec(Gl)
        TTl = TTl * TTT
    end

    ## Create limit cross product matrices

    yyyyd = reshape(GAMM0[:, 1], nobs, nobs) + DD * DD'

    ## cointadd are treated as the first set of variables in XX
    ## coint    are treated as the second set of variables in XX
    ## composition: cointadd - coint - constant - lags
    yyXXd = zeros(S, nobs, p * nobs)
    XXXXd = zeros(S, p * nobs, p * nobs)
    for rr = 1:p
        ## E[yy,x(lag rr)]
        yyXXd[:, nobs * (rr - 1) + 1:nobs * rr] =  reshape(GAMM0[:, rr + 1], nobs, nobs) + DD * DD'

        ## E[x(lag rr),x(lag ll)]
        for ll = rr:p
            yyyydrrll = reshape(GAMM0[:, ll - rr + 1], nobs, nobs) + DD * DD';
            XXXXd[nobs * (rr - 1) + 1:nobs * rr, nobs * (ll - 1) + 1:nobs * ll] =  yyyydrrll
            XXXXd[nobs * (ll - 1) + 1:nobs * ll, nobs * (rr - 1) + 1:nobs * rr] =  yyyydrrll'
        end
    end

    XXyyd = yyXXd'

    if get_VAR
        β = \(XXXXd, XXyyd)
        Σ = yyyyd - XXyyd' * β
        Σ += Σ'  # to correct for machine error
        Σ ./= 2. # and guarantee Σ is symmetric
        return β, Σ
    else
        return yyyyd, XXyyd, XXXXd
    end
end
