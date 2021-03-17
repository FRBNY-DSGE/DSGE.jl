"""
```
measurement(m::AnSchorfheide{T}, TTT::Matrix{T},
            RRR::Matrix{T}, CCC::Vector{T}; reg::Int = 1,
            TTTs::Vector{<: AbstractMatrix{T}} = Matrix{T}[],
            CCCs::Vector{<: AbstractVector{T}} = Vector{T}[],
            information_set::UnitRange = reg:reg,
            memo::Union{ForwardMultipleExpectationsMemo, Nothing} = nothing) where {T <: AbstractFloat}
```

Assign measurement equation

```
y_t = ZZ*s_t + DD + u_t
```

where

```
Var(ϵ_t) = QQ
Var(u_t) = EE
Cov(ϵ_t, u_t) = 0
```
"""
function measurement(m::AnSchorfheide{T},
                     TTT::Matrix{T},
                     RRR::Matrix{T},
                     CCC::Vector{T}; reg::Int = 1,
                     TTTs::Vector{<: AbstractMatrix{T}} = Matrix{T}[],
                     CCCs::Vector{<: AbstractVector{T}} = Vector{T}[],
                     information_set::UnitRange = reg:reg,
                     memo::Union{ForwardMultipleExpectationsMemo, Nothing} = nothing) where {T <: AbstractFloat}

    endo     = m.endogenous_states
    endo_new = m.endogenous_states_augmented
    exo      = m.exogenous_shocks
    obs      = m.observables

    _n_observables = n_observables(m)
    _n_states = n_states_augmented(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    ZZ = zeros(_n_observables, _n_states)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    ## Output growth
    ZZ[obs[:obs_gdp], endo[:y_t]]  = 1.0
    ZZ[obs[:obs_gdp], endo[:y_t1]] = -1.0
    ZZ[obs[:obs_gdp], endo[:z_t]]  = 1.0
    DD[obs[:obs_gdp]]              = m[:γ_Q]

    ## Inflation
    ZZ[obs[:obs_cpi], endo[:π_t]] = 4.0
    DD[obs[:obs_cpi]]             = m[:π_star]

    ## Federal Funds Rate
    ZZ[obs[:obs_nominalrate], endo[:R_t]] = 4.0
    DD[obs[:obs_nominalrate]]             = m[:π_star] + m[:rA] + 4.0*m[:γ_Q]

    # Measurement error
    EE[obs[:obs_gdp], obs[:obs_gdp]]                 = m[:e_y]^2
    EE[obs[:obs_cpi], obs[:obs_cpi]]                 = m[:e_π]^2
    EE[obs[:obs_nominalrate], obs[:obs_nominalrate]] = m[:e_R]^2

    # Variance of innovations
    QQ[exo[:z_sh],exo[:z_sh]]   = (m[:σ_z])^2
    QQ[exo[:g_sh],exo[:g_sh]]   = (m[:σ_g])^2
    QQ[exo[:rm_sh],exo[:rm_sh]] = (m[:σ_R])^2

    return Measurement(ZZ, DD, QQ, EE)
end
