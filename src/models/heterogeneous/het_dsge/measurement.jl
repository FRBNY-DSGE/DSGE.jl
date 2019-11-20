"""
```
measurement(m::SmetsWouters{T}, TTT::Matrix{T}, RRR::Matrix{T},
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
function measurement(m::HetDSGE{T},
                     TTT::Matrix{T},
                     RRR::Matrix{T},
                     CCC::Vector{T}, C_eqn::Vector{T}) where {T<:AbstractFloat}
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
    ZZ[obs[:obs_gdp], first(endo[:y′_t])]       = 1.0
    ZZ[obs[:obs_gdp], first(endo[:y′_t1])] = -1.0
    ZZ[obs[:obs_gdp], first(endo[:z′_t])]       = 1.0
    DD[obs[:obs_gdp]]                   = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    ## Hours growth
    ZZ[obs[:obs_hours], first(endo[:L′_t])] = 1.0
    DD[obs[:obs_hours]]             = m[:Lmean]

    ## Labor Share/real wage growth
    ZZ[obs[:obs_wages], first(endo[:w′_t])]       = 1.0
    ZZ[obs[:obs_wages], first(endo[:w′_t1])] = -1.0
    ZZ[obs[:obs_wages], first(endo[:z′_t])]       = 1.0
    DD[obs[:obs_wages]]                   = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    ## Inflation (GDP Deflator)
    ZZ[obs[:obs_gdpdeflator], first(endo[:π′_t])]  = 1.0
    DD[obs[:obs_gdpdeflator]]              = 100*(m[:π_star]-1)

    ## Nominal interest rate
    ZZ[obs[:obs_nominalrate], first(endo[:R′_t])]       = 1.0
    DD[obs[:obs_nominalrate]]                   = 1 + m[:r] #m[:Rstarn]

    ## Consumption Growth
    ZZ[obs[:obs_consumption], endo_new[:c_t]]       = 1.0
    ZZ[obs[:obs_consumption], endo_new[:c_t1]] = -1.0
    ZZ[obs[:obs_consumption], first(endo[:z′_t])]       = 1.0
    DD[obs[:obs_consumption]]                   = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    ## Investment Growth
    ZZ[obs[:obs_investment], first(endo[:i′_t])]       = 1.0
    ZZ[obs[:obs_investment], endo_new[:i_t1]] = -1.0
    ZZ[obs[:obs_investment], first(endo[:z′_t])]       = 1.0
    DD[obs[:obs_investment]]                   = 100*(exp(m[:γ])-1) #100*(exp(m[:zstar])-1)

    #Measurement error
    EE[obs[:obs_gdp],1] = m[:e_y]^2
    EE[obs[:obs_hours],2] = m[:e_L]^2
    EE[obs[:obs_wages],3] = m[:e_w]^2
    EE[obs[:obs_gdpdeflator],4] = m[:e_π]^2
    EE[obs[:obs_nominalrate],5] = m[:e_R]^2
    EE[obs[:obs_consumption],6] = m[:e_c]^2
    EE[obs[:obs_investment],7] = m[:e_i]^2

    #Variance of innovations
    QQ[exo[:g_sh], exo[:g_sh]]           = m[:σ_g]^2
    QQ[exo[:b_sh], exo[:b_sh]]           = m[:σ_b]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]           = m[:σ_μ]^2
    QQ[exo[:z_sh], exo[:z_sh]]           = m[:σ_z]^2
    QQ[exo[:λ_f_sh], exo[:λ_f_sh]]       = m[:σ_λ_f]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]       = m[:σ_λ_w]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]         = m[:σ_rm]^2

  #=  # These lines set the standard deviations for the anticipated
    # shocks to be equal to the standard deviation for the
    # unanticipated policy shock
    for i = 1:n_anticipated_shocks(m)
        ZZ[obs[Symbol("obs_nominalrate$i")], :] = ZZ[obs[:obs_nominalrate], :]' * (TTT^i)
        DD[obs[Symbol("obs_nominalrate$i")]]    = m[:Rstarn]
        if subspec(m) == "ss1"
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[Symbol("σ_rm")]^2 / n_anticipated_shocks(m)
        else
            QQ[exo[Symbol("rm_shl$i")], exo[Symbol("rm_shl$i")]] = m[Symbol("σ_rm")]^2 / 16
        end
    end
=#
    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end

function construct_consumption_partial(m::AbstractModel, dF2_dRZ::Vector{Float64}, dF2_dWH::Vector{Float64}, dF2_dTT::Vector{Float64})

    c = m[:cstar].value
    μ = m[:μstar].value
    ω = m[:ωstar].value
    ell = m[:lstar].value
    η = m[:η].value
    xgrid::Vector{Float64} = m.grids[:xgrid].points
    xwts::Vector{Float64}  = m.grids[:xgrid].weights
    sgrid::Vector{Float64} = m.grids[:sgrid].points
    swts::Vector{Float64}  = m.grids[:sgrid].weights
    fgrid::Matrix{Float64} = m.grids[:fgrid]
    xswts = kron(swts,xwts)

    ns = get_setting(m, :ns)

    unc = 1 ./ ell .<= repeat(xgrid,ns) .+ η

    denominator = sum(μ .* ω .* c)

    dC_dELL = (μ .* unc .* xswts .*c)' ./ denominator
    dC_dKF = -(xswts .* c)' ./ denominator
    dC_dR = -(xswts .* c)' * dF2_dRZ ./ denominator
    dC_dZ = (xswts .* c)' * dF2_dRZ ./ denominator
    dC_dW = -(xswts .* c)' * dF2_dWH ./ denominator
    dC_dL = -(xswts .* c)' * dF2_dWH ./ denominator
    dC_dT = -(xswts .* c)' * dF2_dTT ./denominator

    return dC_dELL, dC_dKF, dC_dR, dC_dZ, dC_dW, dC_dL, dC_dT
end

function construct_consumption_eqn(m::AbstractModel, TTT_jump::Matrix{Float64}, dF2_dRZ::Vector{Float64}, dF2_dWH::Vector{Float64}, dF2_dTT::Vector{Float64})

    endo_unnorm = m.endogenous_states_unnormalized

    dC_dELL, dC_dKF, dC_dR, dC_dZ, dC_dW, dC_dL, dC_dT = construct_consumption_partial(m, dF2_dRZ, dF2_dWH, dF2_dTT)
    C_eqn = zeros(n_model_states_unnormalized(m))
    C_eqn[endo_unnorm[:l′_t]] = vec(dC_dELL)
    C_eqn[endo_unnorm[:kf′_t]] = vec(dC_dKF)
    C_eqn[first(endo_unnorm[:R′_t])] = dC_dR
    C_eqn[first(endo_unnorm[:z′_t])] = dC_dZ
    C_eqn[first(endo_unnorm[:w′_t])] = dC_dW
    C_eqn[first(endo_unnorm[:L′_t])] = dC_dL
    C_eqn[first(endo_unnorm[:t′_t])] = dC_dT

    Qx, Qy, _, _ = compose_normalization_matrices(m)

    gx2 = Qy'*TTT_jump*Qx
    n_backward_looking_states_unnorm = n_backward_looking_states_unnormalized(m)

    C = C_eqn'*[Matrix{Float64}(I, n_backward_looking_states_unnorm, n_backward_looking_states_unnorm); gx2]*Qx'
    return vec(C)
end
