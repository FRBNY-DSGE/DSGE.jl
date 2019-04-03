"""
```
measurement{T<:AbstractFloat}(m::RealBond{T}, TTT::Matrix{T}, RRR::Matrix{T},
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
function measurement(m::RealBond{T}, TTT::Matrix{T},
                     TTT_jump::Matrix{T},
                     RRR::Matrix{T}, CCC::Vector{T}) where T<:AbstractFloat
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

    # This is just a snippet of code that constructs linearized GDP
    # It should be incorporated into a larger file that simulates, filters etc.
    # This is based on real_bond_hank.jl

    dGDP_dMU, GDP, dGDP_dELL, dGDP_dRR, dGDP_dWW, dGDP_dTT = construct_GDPfn_realbond(nx, ns, xgrid_total, sgrid, xwts, swts,
                                                                                      γ, ν, abar, R, aborrow, μ,
                                                                                      c, η, ell, χss)

    GDPfn = zeros(1, 2*nx*ns+7) # GDP as function of un-normalized MU Z MON ELL RR II WW PI TT
    # note: GDP is only a function of contemporaneous variables
    # we are using the indices (MUP,ZP, etc.) corresponding to date t+1 variables
    # simply because these indices happen to work here also
    # this does not mean that GDP is a function of date t+1 variables
    GDPfn[1, endo[:μ′_t]]  = dGDP_dMU
    GDPfn[1, endo[:z′_t]]  .= GDP
    GDPfn[1, endo[:l′_t]]  = dGDP_dELL
    GDPfn[1, endo[:R′_t]]  .= dGDP_dRR
    GDPfn[1, endo[:w′_t]]  .= dGDP_dWW
    GDPfn[1, endo[:t′_t]]  .= dGDP_dTT

    ########################################
    Qx, Qy, _, _ = compose_normalization_matrices(m)
    gx2  = Qy'*TTT_jump*Qx

    # now we need to create GDP as a function of the normalized states
    ZZ_states = (1/GDP)*GDPfn*[eye(nx*ns+2); gx2]*Qx' # this is for log GDP
                                                      # to use the level of gdp, remove (1/GDP)

    ZZ = Matrix{Float64}(undef, _n_observables, _n_model_states)
    ZZ[1:_n_states] = ZZ_states
    ZZ[_n_states+1:end] = zeros(_n_jumps)

    # Measurement error
    EE[obs[:obs_gdp], obs[:obs_gdp]] = m[:e_y]

    # Variance of innovations
    QQ[exo[:z_sh], exo[:z_sh]] = m[:σ_z]^2
    QQ[exo[:mon_sh], exo[:mon_sh]] = m[:σ_mon]^2

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end

function construct_GDPfn_realbond(nx::Int, ns::Int,
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
    GDP       = 0.0

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
            GDP         += swts[is]*xwts[ix]*sgrid[is]*η[i]*μ[i] # note, this is also dGDP_dZ. It is already created by real_bond_hank.jl
        end
    end
    return dGDP_dMU, GDP, dGDP_dELL, dGDP_dRR, dGDP_dWW, dGDP_dTT
end
