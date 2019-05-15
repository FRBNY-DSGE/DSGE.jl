"""
```
measurement(m::RepDSGEGovDebt{T}, TTT::Matrix{T}, RRR::Matrix{T},
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
function measurement(m::RepDSGEGovDebt{T},
                     TTT::Matrix{T},
                     RRR::Matrix{T},
                     CCC::Vector{T}) where {T<:AbstractFloat}
    endo      = m.endogenous_states
    endo_new  = m.endogenous_states_augmented
    exo       = m.exogenous_shocks
    obs       = m.observables

    _n_model_states = get_setting(m, :n_model_states_augmented)
    _n_states = n_backward_looking_states(m)
    _n_jumps = n_jumps(m)

    _n_observables = n_observables(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    ZZ = zeros(_n_observables, _n_model_states)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    ## Output growth - Quarterly!
    ZZ[obs[:obs_gdp], endo[:y′_t]]  = 1.0
    ZZ[obs[:obs_gdp], endo[:y′_t1]] = -1.0
    ZZ[obs[:obs_gdp], endo[:z′_t]]  = 1.0
    DD[obs[:obs_gdp]]               = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    ## Hours growth
    ZZ[obs[:obs_hours], endo[:L′_t]] = 1.0
    DD[obs[:obs_hours]]              = m[:Lmean]

    ## Labor Share/real wage growth
    ZZ[obs[:obs_wages], endo[:w′_t]]       = 1.0
    ZZ[obs[:obs_wages], endo[:w′_t1]]      = -1.0
    ZZ[obs[:obs_wages], endo[:z′_t]]       = 1.0
    DD[obs[:obs_wages]]                    = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    ## Inflation (GDP Deflator)
    ZZ[obs[:obs_gdpdeflator], endo[:π′_t]] = 1.0
    DD[obs[:obs_gdpdeflator]]              = 100*(m[:π_star]-1)

    ## Nominal interest rate
    ZZ[obs[:obs_nominalrate], endo[:R′_t]] = 1.0
    DD[obs[:obs_nominalrate]]              = 1 + m[:r] #m[:Rstarn]

    ## Consumption Growth
    ZZ[obs[:obs_consumption], endo[:l′_t]]        = -1.0
    ZZ[obs[:obs_consumption], endo_new[:c_t1]]    = -1.0
    ZZ[obs[:obs_consumption], endo[:z′_t]]        = 1.0
    DD[obs[:obs_consumption]]                     = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    ## Investment Growth
    ZZ[obs[:obs_investment], endo[:i′_t]]       = 1.0
    ZZ[obs[:obs_investment], endo_new[:i_t1]]   = -1.0
    ZZ[obs[:obs_investment], endo[:z′_t]]       = 1.0
    DD[obs[:obs_investment]]                    = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    # Measurement error
    EE[obs[:obs_gdp],1]         = m[:e_y]^2
    EE[obs[:obs_hours],2]       = m[:e_L]^2
    EE[obs[:obs_wages],3]       = m[:e_w]^2
    EE[obs[:obs_gdpdeflator],4] = m[:e_π]^2
    EE[obs[:obs_nominalrate],5] = m[:e_R]^2
    EE[obs[:obs_consumption],6] = m[:e_c]^2
    EE[obs[:obs_investment],7]  = m[:e_i]^2

    #Variance of innovations
    QQ[exo[:g_sh], exo[:g_sh]]           = m[:σ_g]^2
    QQ[exo[:b_sh], exo[:b_sh]]           = m[:σ_b]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]           = m[:σ_μ]^2
    QQ[exo[:z_sh], exo[:z_sh]]           = m[:σ_z]^2
    QQ[exo[:λ_f_sh], exo[:λ_f_sh]]       = m[:σ_λ_f]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]       = m[:σ_λ_w]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]         = m[:σ_rm]^2

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end
