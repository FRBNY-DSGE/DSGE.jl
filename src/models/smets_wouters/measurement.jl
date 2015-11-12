# Assign measurement equation : X_t = ZZ*S_t + DD + u_t
# where u_t = eta_t+MM* eps_t with var(eta_t) = EE
# where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
function measurement{T<:AbstractFloat}(m::SmetsWouters{T},
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
        _num_observables = num_observables(m)
        _num_states = num_states_augmented(m)
        _num_shocks_exogenous = num_shocks_exogenous(m)
        endo_addl = m.endogenous_states_postgensys
    else
        _num_observables = num_observables(m) - num_anticipated_shocks(m)
        _num_states = num_states_augmented(m) - num_anticipated_shocks(m)
        _num_shocks_exogenous = num_shocks_exogenous(m) - num_anticipated_shocks(m)
        endo_addl = Dict(
            [(key,m.endogenous_states_postgensys[key] - num_anticipated_shocks(m)) for key in keys(m.endogenous_states_postgensys)])
    end

    ZZ = zeros(_num_observables, _num_states)
    DD = zeros(_num_observables, 1)
    MM = zeros(_num_observables, _num_shocks_exogenous)
    EE = zeros(_num_observables, _num_observables)
    QQ = zeros(_num_shocks_exogenous, _num_shocks_exogenous)

    ## Output growth - Quarterly!
    ZZ[obs[:g_y], endo[:y_t]]       = 1.0
    ZZ[obs[:g_y], endo_addl[:y_t1]] = -1.0
    ZZ[obs[:g_y], endo[:z_t]]       = 1.0
    DD[obs[:g_y]]                   = 100*(exp(m[:zstar])-1)

    ## Hours growth
    ZZ[obs[:g_hours], endo[:L_t]] = 1.0
    DD[obs[:g_hours]]             = m[:Lmean]

    ## Labor Share/real wage growth
    ZZ[obs[:g_w], endo[:w_t]]       = 1.0
    ZZ[obs[:g_w], endo_addl[:w_t1]] = -1.0
    ZZ[obs[:g_w], endo[:z_t]]       = 1.0
    DD[obs[:g_w]]                   = 100*(exp(m[:zstar])-1)

    ## Inflation (GDP Deflator)
    ZZ[obs[:π_gdpdef], endo[:π_t]]  = 1.0
    DD[obs[:π_gdpdef]]              = 100*(m[:π_star]-1)

    ## Nominal interest rate
    ZZ[obs[:R_n], endo[:R_t]]       = 1.0
    DD[obs[:R_n]]                   = m[:Rstarn]

    ## Consumption Growth
    ZZ[obs[:g_c], endo[:c_t]]       = 1.0
    ZZ[obs[:g_c], endo_addl[:c_t1]] = -1.0
    ZZ[obs[:g_c], endo[:z_t]]       = 1.0
    DD[obs[:g_c]]                   = 100*(exp(m[:zstar])-1)

    ## Investment Growth
    ZZ[obs[:g_i], endo[:i_t]]       = 1.0
    ZZ[obs[:g_i], endo_addl[:i_t1]] = -1.0
    ZZ[obs[:g_i], endo[:z_t]]       = 1.0
    DD[obs[:g_i]]                   = 100*(exp(m[:zstar])-1)

    QQ[exo[:g_sh], exo[:g_sh]]           = m[:σ_g]^2
    QQ[exo[:b_sh], exo[:b_sh]]           = m[:σ_b]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]           = m[:σ_μ]^2
    QQ[exo[:z_sh], exo[:z_sh]]           = m[:σ_z]^2
    QQ[exo[:λ_f_sh], exo[:λ_f_sh]]       = m[:σ_λ_f]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]       = m[:σ_λ_w]^2
    QQ[exo[:rm_sh], exo[:rm_sh]]         = m[:σ_rm]^2

    # These lines set the standard deviations for the anticipated
    # shocks to be equal to the standard deviation for the
    # unanticipated policy shock
    if shocks
        for i = 1:num_anticipated_shocks(m)
            QQ[exo[symbol("rm_shl$i")], exo[symbol("rm_shl$i")]] = m[symbol("σ_rm$i")]^2
        end
    end

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC != 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    HH = EE + MM*QQ*MM'
    VV = QQ*MM'
    VVall = [[RRR*QQ*RRR' RRR*VV];
             [VV'*RRR'    HH]]

    return Measurement(ZZ, DD, QQ, EE, MM, VVall)
end
