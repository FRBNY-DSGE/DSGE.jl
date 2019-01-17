import DSGE: Measurement, get_setting, n_observables, n_shocks_exogenous

"""
```
measurement(m::KrusellSmith{T}, TTT::Matrix{T}, RRR::Matrix{T},
            CCC::Vector{T}) where {T<:AbstractFloat}
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
function measurement(m::KrusellSmith{T},
                     TTT::Matrix{T}, TTT_jump::Matrix{T},
                     RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
    endo      = m.endogenous_states
    exo       = m.exogenous_shocks
    obs       = m.observables

    _n_model_states = n_model_states(m)

    _n_observables = n_observables(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    ZZ = zeros(_n_observables, _n_model_states)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    # Output in log levels
    ZZ[obs[:obs_gdp], endo[:z′_t]]   = (m[:Kstar]^m[:α])*(m[:Lstar]^(1.0-m[:α]))
    ZZ[obs[:obs_gdp], endo[:K′_t]]   = (m[:α]/m[:Kstar])*(m[:Kstar]^m[:α])*(m[:Lstar]^(1.0-m[:α]))

    # Measurement error
    EE[obs[:obs_gdp], obs[:obs_gdp]] = m[:e_y]

    # Variance of innovations
    QQ[exo[:z_sh], exo[:z_sh]] = m[:σ_z]^2

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC != 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end
