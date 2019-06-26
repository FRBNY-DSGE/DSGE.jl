"""
```
measurement{T<:AbstractFloat}(m::BondLabor{T}, TTT::Matrix{T}, RRR::Matrix{T},
                              CCC::Vector{T})
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
function measurement(m::BondLabor{T}, TTT::Matrix{T},
                     TTT_jump::Matrix{T},
                     RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
    endo      = m.endogenous_states_unnormalized
    exo       = m.exogenous_shocks
    obs       = m.observables

    _n_model_states = n_model_states(m)
    _n_states       = n_backward_looking_states(m)
    _n_jumps        = n_jumps(m)

    _n_observables = n_observables(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    ZZ = zeros(_n_observables, _n_model_states)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    nx = get_setting(m, :nx)
    ns = get_setting(m, :ns)

    ν   = m[:ν]
    γ   = m[:γ]
    abar = m[:abar]
    R    = m[:R]

    l   = m[:lstar].value
    μ   = m[:μstar].value
    c   = m[:cstar].value
    η   = m[:ηstar].value
    χss = m[:χstar].value

    xgrid = m.grids[:xgrid].points
    xgrid_total = m.grids[:xgrid_total]
    xwts  = m.grids[:xgrid].weights
    sgrid = m.grids[:sgrid].points
    sgrid_total = m.grids[:sgrid_total]
    swts  = m.grids[:sgrid].weights
    weights_total = m.grids[:weights_total]

    aborrow = abar/R
    chipW = (1+ν)./(ν./(aborrow + χss - xgrid_total) + γ./χss)
    chipR = (1/(R*R))*(ν*abar./(aborrow + χss - xgrid_total))./(ν./(aborrow + χss - xgrid_total) + γ./χss)

    GDP  = 0.0
    GDPZ = 0.0
    GDPR = 0.0

    GDPfn = zeros(1,2*ns*nx +2) # GDP as function of un-normalized MU Z l R
    for ix =1:nx
        for is=1:ns
            i = ix + nx*(is-1)
            GDP += μ[i]*swts[is]*xwts[ix]*η[i]
            GDPZ += sgrid[is]*μ[i]*swts[is]*xwts[ix]*(η[i]/ν - ((γ*η[i])/(ν*c[i]))*chipW[i]*((l[i]^(-1/γ))>χss[i]) )
            GDPR += -(sgrid[is]*μ[i]*γ*η[i]/(ν*c[i]))*swts[is]*xwts[ix]*chipR[i]*((l[i]^(-1/γ))>χss[i])
        end
    end

    # Output in log levels
    GDPfn[1, endo[:μ′_t]] = η'.*sgrid_total'.*weights_total'
    GDPfn[1, endo[:z′_t]] = GDPZ + GDP
    GDPfn[1, endo[:l′_t]] = (1/γ)*μ'.*(sgrid_total'*(γ/ν).*(η./c)').*weights_total'.*( l.^(-1-1/γ) )'.*((l.^(-1/γ)).<=χss)'
    GDPfn[1, endo[:R′_t]] = GDPR

    Qx, Qy, _, _ = compose_normalization_matrices(m)
    gx2  = Qy'*TTT_jump*Qx

    # now we need to create GDP as a function of the normalized states
    ZZ_states = (1/GDP)*GDPfn*[eye(nx*ns+1); gx2]*Qx' # this is for log GDP
                                                      # to use the level of gdp, remove (1/GDP)
    ZZ = Matrix{Float64}(_n_observables, _n_model_states)
    ZZ[1:_n_states] = ZZ_states
    ZZ[_n_states+1:end] = zeros(_n_jumps)

    # Measurement error
    EE[obs[:obs_gdp], obs[:obs_gdp]] = m[:e_y]

    # Variance of innovations
    QQ[exo[:z_sh], exo[:z_sh]] = m[:σ_z]^2

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end
