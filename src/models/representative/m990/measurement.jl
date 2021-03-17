"""
```
measurement(m::Model990{T}, TTT::Matrix{T}, RRR::Matrix{T},
            CCC::Vector{T}; reg::Int = 1,
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
function measurement(m::Model990{T}, TTT::Matrix{T}, RRR::Matrix{T},
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

    ## Output growth - Quarterly!
    ZZ[obs[:obs_gdp], endo[:y_t]]      = 1.0
    ZZ[obs[:obs_gdp], endo_new[:y_t1]] = -1.0
    ZZ[obs[:obs_gdp], endo[:z_t]]      = 1.0
    DD[obs[:obs_gdp]]                  = 100*(exp(m[:z_star])-1)

    ## Hours growth
    ZZ[obs[:obs_hours], endo[:L_t]] = 1.0
    DD[obs[:obs_hours]]             = m[:Lmean]

    ## Labor Share/real wage growth
    ZZ[obs[:obs_wages], endo[:w_t]]      = 1.0
    ZZ[obs[:obs_wages], endo_new[:w_t1]] = -1.0
    ZZ[obs[:obs_wages], endo[:z_t]]      = 1.0
    DD[obs[:obs_wages]]                  = 100*(exp(m[:z_star])-1)

    ## Inflation (GDP Deflator)
    ZZ[obs[:obs_gdpdeflator], endo[:π_t]]            = m[:Γ_gdpdef]
    ZZ[obs[:obs_gdpdeflator], endo_new[:e_gdpdef_t]] = 1.0
    DD[obs[:obs_gdpdeflator]]                        = 100*(m[:π_star]-1) + m[:δ_gdpdef]

    ## Inflation (Core PCE)
    ZZ[obs[:obs_corepce], endo[:π_t]]             = 1.0
    ZZ[obs[:obs_corepce], endo_new[:e_corepce_t]] = 1.0
    DD[obs[:obs_corepce]]                         = 100*(m[:π_star]-1)

    ## Nominal interest rate
    ZZ[obs[:obs_nominalrate], endo[:R_t]] = 1.0
    DD[obs[:obs_nominalrate]]             = m[:Rstarn]

    ## Consumption Growth
    ZZ[obs[:obs_consumption], endo[:c_t]]      = 1.0
    ZZ[obs[:obs_consumption], endo_new[:c_t1]] = -1.0
    ZZ[obs[:obs_consumption], endo[:z_t]]      = 1.0
    DD[obs[:obs_consumption]]                  = 100*(exp(m[:z_star])-1)

    ## Investment Growth
    ZZ[obs[:obs_investment], endo[:i_t]]      = 1.0
    ZZ[obs[:obs_investment], endo_new[:i_t1]] = -1.0
    ZZ[obs[:obs_investment], endo[:z_t]]      = 1.0
    DD[obs[:obs_investment]]                  = 100*(exp(m[:z_star])-1)

    ## Spreads
    ZZ[obs[:obs_spread], endo[:ERtil_k_t]] = 1.0
    ZZ[obs[:obs_spread], endo[:R_t]]       = -1.0
    DD[obs[:obs_spread]]                   = 100*log(m[:spr])

    ## 10 yrs infl exp
    TTT10                          = (1/40)*((Matrix{Float64}(I, size(TTT, 1), size(TTT, 1))
                                              - TTT)\(TTT - TTT^41))
    ZZ[obs[:obs_longinflation], :] = TTT10[endo[:π_t], :]
    DD[obs[:obs_longinflation]]    = 100*(m[:π_star]-1)

    ## Long Rate
    ZZ[obs[:obs_longrate], :]               = ZZ[6, :]' * TTT10
    ZZ[obs[:obs_longrate], endo_new[:lr_t]] = 1.0
    DD[obs[:obs_longrate]]                  = m[:Rstarn]

    ## TFP
    ZZ[obs[:obs_tfp], endo[:z_t]]       = (1-m[:α])*m[:Iendoα] + 1*(1-m[:Iendoα])
    ZZ[obs[:obs_tfp], endo_new[:tfp_t]] = 1.0
    ZZ[obs[:obs_tfp], endo[:u_t]]       = m[:α]/( (1-m[:α])*(1-m[:Iendoα]) + 1*m[:Iendoα] )
    ZZ[obs[:obs_tfp], endo_new[:u_t1]]  = -(m[:α]/( (1-m[:α])*(1-m[:Iendoα]) + 1*m[:Iendoα]) )

    QQ[exo[:g_sh], exo[:g_sh]]             = m[:σ_g]^2
    QQ[exo[:b_sh], exo[:b_sh]]             = m[:σ_b]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]             = m[:σ_μ]^2
    QQ[exo[:z_sh], exo[:z_sh]]             = m[:σ_z]^2
    QQ[exo[:λ_f_sh], exo[:λ_f_sh]]         = m[:σ_λ_f]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]         = m[:σ_λ_w]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]           = m[:σ_r_m]^2
    QQ[exo[:σ_ω_sh], exo[:σ_ω_sh]]         = m[:σ_σ_ω]^2
    QQ[exo[:μ_e_sh], exo[:μ_e_sh]]         = m[:σ_μ_e]^2
    QQ[exo[:γ_sh], exo[:γ_sh]]             = m[:σ_γ]^2
    QQ[exo[:π_star_sh], exo[:π_star_sh]]   = m[:σ_π_star]^2
    QQ[exo[:lr_sh], exo[:lr_sh]]           = m[:σ_lr]^2
    QQ[exo[:zp_sh], exo[:zp_sh]]           = m[:σ_z_p]^2
    QQ[exo[:tfp_sh], exo[:tfp_sh]]         = m[:σ_tfp]^2
    QQ[exo[:gdpdef_sh], exo[:gdpdef_sh]]   = m[:σ_gdpdef]^2
    QQ[exo[:corepce_sh], exo[:corepce_sh]] = m[:σ_corepce]^2

    # These lines set the standard deviations for the anticipated shocks
    for i = 1:n_mon_anticipated_shocks(m)
        ZZ[obs[Symbol("obs_nominalrate$i")], :] = ZZ[obs[:obs_nominalrate], :]' * (TTT^i)
        DD[obs[Symbol("obs_nominalrate$i")]]    = m[:Rstarn]
        if subspec(m) == "ss6"
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[:σ_r_m]^2 / n_mon_anticipated_shocks(m)
        else
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[Symbol("σ_r_m$i")]^2
        end
    end

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end
