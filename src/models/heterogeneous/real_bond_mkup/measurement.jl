"""
```
measurement{T<:AbstractFloat}(m::RealBondMkup{T}, TTT::Matrix{T}, RRR::Matrix{T},
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
function measurement(m::RealBondMkup{T}, TTT::Matrix{T},
                     TTT_jump::Matrix{T},
                     RRR::Matrix{T}, CCC::Vector{T}, GDPeqn::AbstractArray{T}) where {T<:AbstractFloat}

    #@info "enter"
    endo        = m.endogenous_states
    endo_unnorm = m.endogenous_states_unnormalized
    endo_new    = m.model_states_augmented
    exo         = m.exogenous_shocks
    obs         = m.observables

    _n_model_states = n_model_states(m)
    _n_model_states_aug = _n_model_states + length(endo_new)
    _n_states       = n_backward_looking_states(m)
    _n_jumps        = n_jumps(m)

    _n_observables = n_observables(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    # Load in parameters, steady-state parameters, and grids
    R::Float64     = m[:R].value

    ZZ = zeros(_n_observables, _n_model_states_aug)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    # GDP
    ZZ[obs[:obs_gdp], 1:_n_states]     = GDPeqn
    ZZ[obs[:obs_gdp], _n_states+1:_n_states + _n_jumps] = zeros(_n_jumps)
    #@show size(ZZ[obs[:obs_gdp], _n_states + _n_jumps + 1:end])
    ZZ[obs[:obs_gdp], first(_n_states + _n_jumps + 1:end)] = -1

    # Inflation
    ZZ[obs[:obs_corepce], first(endo[:π′_t])] = 1.

    # Nominal FFR
    ZZ[obs[:obs_nominalrate], first(endo[:i′_t])] = 1.0 / R

    # Measurement error
    EE[obs[:obs_gdp], obs[:obs_gdp]] = m[:e_y]

    # Variance of innovations
    QQ[exo[:z_sh], exo[:z_sh]] = m[:σ_z]^2
    QQ[exo[:mon_sh], exo[:mon_sh]] = m[:σ_mon]^2
    QQ[exo[:mkp_sh], exo[:mkp_sh]] = m[:σ_mkp]^2

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0)
        @show CCC[CCC.!=0]
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end

function construct_GDPfn_realbondmkp(nx::Int, ns::Int,
                                  xgrid_total::Vector{Float64}, sgrid::Vector{Float64},
                                  xwts::Vector{Float64}, swts::Vector{Float64},
                                  γ::Float64, ν::Float64, abar::Float64,
                                  R::Float64, aborrow::Float64,
                                  μ::Vector{Float64}, c::Vector{Float64}, η::Vector{Float64},
                                  ell::Vector{Float64}, χss::Vector{Float64})
    dGDP_dMU  = zeros(1,nx*ns)
    dGDP_dELL = zeros(1,nx*ns)
    dGDP_dWW  = 0.0
    dGDP_dRR  = 0.0
    dGDP_dTT  = 0.0
    dGDP_dZ   = 0.0

    unc       = zeros(nx*ns)
    chipW, chipR, chipX = construct_chip_realbond(xgrid_total, γ, ν, aborrow,
                                                  abar, R, χss)
    for ix =1:nx
        for is=1:ns
            i = ix + nx*(is-1)
            unc[i] = ((ell[i]^(-1/γ))<=χss[i]) # =1 if unconstrained
            dGDP_dMU[i]  = xwts[ix]*swts[is]*sgrid[is]*η[i]
            dGDP_dELL[i] = xwts[ix]*swts[is]*μ[ix]*(sgrid[is]*γ*η[i]/(ν*c[i]))*(1/γ)*unc[i]*ell[i]^(-(1/γ)-1)
            dGDP_dWW    += swts[is]*xwts[ix]*μ[i]*((sgrid[is]*η[i])/ν - (sgrid[is]*γ*η[i]/(ν*c[i]))*(1-unc[i])*chipW[i])
            dGDP_dRR    += -swts[is]*xwts[ix]*μ[i]*(sgrid[is]*γ*η[i]/(ν*c[i]))*(1-unc[i])*chipR[i]
            dGDP_dTT    += -swts[is]*xwts[ix]*μ[i]*(sgrid[is]*γ*η[i]/(ν*c[i]))*(1-unc[i])*chipX[i]
            dGDP_dZ     += swts[is]*xwts[ix]*sgrid[is]*η[i]*μ[i]
        end
    end
    return dGDP_dMU, dGDP_dZ, dGDP_dELL, dGDP_dRR, dGDP_dWW, dGDP_dTT
end

function update_measurement_covariance_matrices!(m::RealBondMkup, system::System{T}) where {T<:AbstractFloat}
    exo       = m.exogenous_shocks

    # Variance of innovations
    system[:QQ][exo[:z_sh], exo[:z_sh]] = m[:σ_z]^2
    system[:QQ][exo[:mon_sh], exo[:mon_sh]] = m[:σ_mon]^2
    system[:QQ][exo[:mkp_sh], exo[:mkp_sh]] = m[:σ_mkp]^2
end

function construct_GDPeqn(m::RealBondMkup, TTT_jump::Matrix{T}) where {T<:AbstractFloat}
    endo_unnorm = m.endogenous_states_unnormalized

    # Load in parameters, steady-state parameters, and grids
    γ::Float64     = m[:γ].value
    ν::Float64     = m[:ν].value
    R::Float64     = m[:R].value
    abar::Float64  = m[:abar].value
    aborrow        = abar/R

    ell::Vector{Float64}  = m[:lstar].value
    c::Vector{Float64}    = m[:cstar].value
    μ::Vector{Float64}    = m[:μstar].value
    η::Vector{Float64}    = m[:ηstar].value
    χss::Vector{Float64}  = m[:χstar].value

    xwts::Vector{Float64}  = m.grids[:xgrid].weights
    sgrid::Vector{Float64} = m.grids[:sgrid].points
    swts::Vector{Float64}  = m.grids[:sgrid].weights
    xgrid_total::Vector{Float64}  = m.grids[:xgrid_total]

    nx::Int = get_setting(m, :nx)
    ns::Int = get_setting(m, :ns)

    # Construct GDP
    dGDP_dMU, dGDP_dZ, dGDP_dELL, dGDP_dRR, dGDP_dWW, dGDP_dTT =
        construct_GDPfn_realbondmkp(nx, ns, xgrid_total, sgrid, xwts, swts,
                                 γ, ν, abar, R, aborrow, μ, c, η, ell, χss)

    GDPfn = zeros(n_model_states_unnormalized(m))
    # GDP as function of un-normalized MU Z MON ELL RR II WW PI TT
    # note: GDP is only a function of contemporaneous variables
    # we are using the indices (MUP,ZP, etc.) corresponding to date t+1 variables
    # simply because these indices happen to work here also
    # this does not mean that GDP is a function of date t+1 variables
    GDPfn[endo_unnorm[:μ′_t]]  = vec(dGDP_dMU)
    GDPfn[first(endo_unnorm[:z′_t])]  = dGDP_dZ      #first because the index is a range of single number
    GDPfn[endo_unnorm[:l′_t]]  = vec(dGDP_dELL)
    GDPfn[first(endo_unnorm[:R′_t])]  = dGDP_dRR
    GDPfn[first(endo_unnorm[:w′_t])]  = dGDP_dWW
    GDPfn[first(endo_unnorm[:t′_t])]  = dGDP_dTT

    ########################################

    Qx, Qy, _, _ = compose_normalization_matrices(m)

    gx2  = Qy'*TTT_jump*Qx

    # now we need to create GDP as a function of the normalized states
    n_backward_looking_states_unnorm = n_backward_looking_states_unnormalized(m)
    GDPeqn = GDPfn'*[Matrix{Float64}(I, n_backward_looking_states_unnorm, n_backward_looking_states_unnorm); gx2]*Qx'
             # this is for level of GDP
             # to use the log of gdp, front multiply by (1/dGDP_dZ)
    return GDPeqn
end
