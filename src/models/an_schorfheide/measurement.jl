"""
```
measurement{T<:AbstractFloat}(m::AnSchorfheide{T}, TTT::Matrix{T}, RRR::Matrix{T},
                              CCC::Matrix{T}; shocks::Bool = true)
```

Assign measurement equation
```
X_t = ZZ*S_t + DD + u_t
```
where
```
u_t = eta_t + MM*eps_t
var(eta_t) = EE
var(u_t) = HH = EE + MM*QQ*MM'
cov(eps_t,u_t) = VV = QQ*MM'
```
"""
function measurement{T<:AbstractFloat}(m::AnSchorfheide{T},
                                       TTT::Matrix{T},
                                       RRR::Matrix{T},
                                       CCC::Matrix{T};
                                       shocks::Bool = true)
endo = m.endogenous_states
exo  = m.exogenous_shocks
obs  = m.observables

# If shocks = true, then return measurement equation matrices with rows and columns for
# anticipated policy shocks
if shocks
    _n_observables = n_observables(m)
    _n_states = n_states_augmented(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)
    endo_new = m.endogenous_states_augmented
else
    _n_observables = n_observables(m) - n_anticipated_shocks(m)
    _n_states = n_states_augmented(m) - n_anticipated_shocks(m)
    _n_shocks_exogenous = n_shocks_exogenous(m) - n_anticipated_shocks(m)
    endo_new = Dict(
        [(key,m.endogenous_states_augmented[key] - n_anticipated_shocks(m)) for key in keys(m.endogenous_states_augmented)])
end

    ZZ = zeros(_n_observables, _n_states)
    DD = zeros(_n_observables, 1)
    MM = zeros(_n_observables, _n_shocks_exogenous)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)


    ## Output growth
    ZZ[obs[:obs_gdp], endo[:y_t]]       = 1.0
    ZZ[obs[:obs_gdp], endo[:y_t1]]      = -1.0
    ZZ[obs[:obs_gdp], endo[:z_t]]       = 1.0
    DD[obs[:obs_gdp]]                   = m[:γ_Q]

    ## Inflation
    ZZ[obs[:obs_infl], endo[:π_t]]   = 4.0
    DD[obs[:obs_infl]]               = m[:π_star]

    ## Federal Funds Rate
    ZZ[obs[:obs_ffr], endo[:R_t]]       = 4.0
    DD[obs[:obs_ffr]]                   = m[:π_star] + m[:rA] + 4.0*m[:γ_Q]

    # Measurement error
    EE[obs[:obs_gdp], endo[:y_t]]     = m[:e_y]^2
    EE[obs[:obs_infl], endo[:π_t]] = m[:e_π]^2
    EE[obs[:obs_ffr], endo[:R_t]]     = m[:e_R]^2

    # Variance of innovations
    QQ[exo[:z_sh],exo[:z_sh]] = (m[:σ_z])^2
    QQ[exo[:g_sh],exo[:g_sh]] = (m[:σ_g])^2
    QQ[exo[:R_sh],exo[:R_sh]] = (m[:σ_R])^2

    HH = EE + MM*QQ*MM'
    VV    = QQ*MM'
    VVall = [[RRR*QQ*RRR' RRR*VV];
             [VV'*RRR'    HH]]

    return Measurement(ZZ, DD, QQ, EE, MM, VVall)
end
