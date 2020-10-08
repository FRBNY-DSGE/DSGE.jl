"""
```
taylor_replace_eq_entries(m::AbstractDSGEModel,
                                         Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                                         C::Vector{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64})taylor_replace_eq_entries()
```

implements a Taylor policy rule with monetary policy shocks and
"""

function taylor_replace_eq_entries(m::AbstractDSGEModel,
                                         Γ0::Matrix{Float64}, Γ1::Matrix{Float64},
                                         C::Vector{Float64}, Ψ::Matrix{Float64}, Π::Matrix{Float64})
    # add law of motion for pgap
    eq        = m.equilibrium_conditions
    endo      = m.endogenous_states
#=
    ait_Thalf = haskey(get_settings(m), :ait_Thalf) ? get_setting(m, :ait_Thalf) : 10.
    gdp_Thalf = haskey(get_settings(m), :gdp_Thalf) ? get_setting(m, :gdp_Thalf) : 10.
    ρ_pgap    = exp(log(0.5) / ait_Thalf)
    ρ_ygap    = exp(log(0.5) / gdp_Thalf)
    ρ_smooth  = haskey(get_settings(m), :flexible_ait_ρ_smooth) ? get_setting(m, :flexible_ait_ρ_smooth) : 0.656 # m[:ρ]
    φ_π       = haskey(get_settings(m), :flexible_ait_φ_π) ? get_setting(m, :flexible_ait_φ_π) : 11.13
    φ_y       = haskey(get_settings(m), :flexible_ait_φ_y) ? get_setting(m, :flexible_ait_φ_y) : 11.13

    # This assumes that the inflation target is the model's steady state
    Γ0[eq[:eq_pgap], endo[:pgap_t]] = 0.
    Γ0[eq[:eq_pgap], endo[:π_t]]    = 0.
    Γ1[eq[:eq_pgap], endo[:pgap_t]] = 0.

    # Add in the GDP targeting rule
    Γ0[eq[:eq_ygap], endo[:ygap_t]] = 0.
    # Γ0[eq[:eq_ygap], endo[:π_t]]    = -1. # already accounted for by AIT
    Γ1[eq[:eq_ygap], endo[:ygap_t]] = 0. # 1.

    Γ0[eq[:eq_ygap], endo[:y_t]]    = 0.
    Γ0[eq[:eq_ygap], endo[:z_t]]    = 0.
    Γ1[eq[:eq_ygap], endo[:y_t]]    = 0.
=#
    # Zero out old policy rule
    Γ0[eq[:eq_mp], :]              .= 0.
    Γ1[eq[:eq_mp], :]              .= 0.
    Ψ[eq[:eq_mp], :]               .= 0.

    # Replace monetary policy rule
    Γ0[eq[:eq_mp], endo[:R_t]]      = 1.
    Γ1[eq[:eq_mp], endo[:R_t]]      = m[:ρ]
    C[ eq[:eq_mp]]                  = 0.

    Γ0[eq[:eq_mp], endo[:π_t]]      = -(1-m[:ρ])*m[:ψ1]
    Γ0[eq[:eq_mp], endo[:π_star_t]] = (1-m[:ρ])*m[:ψ1]
    Γ0[eq[:eq_mp], endo[:y_t]]      = -(1-m[:ρ])*m[:ψ2] - m[:ψ3]
    Γ0[eq[:eq_mp], endo[:y_f_t]]    = (1-m[:ρ])*m[:ψ2] + m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_t]]      = -m[:ψ3]
    Γ1[eq[:eq_mp], endo[:y_f_t]]    = m[:ψ3]
    # Γ0[eq[:eq_mp], endo[:rm_t]]      = -1.

    # Γ0[eq[:eq_mp], endo[:pgap_t]]   = -φ_π * (1. - ρ_pgap) * (1. - ρ_smooth) # (1 - m[:ρ]) # This is the AIT part
    # Γ0[eq[:eq_mp], endo[:ygap_t]]   = -φ_y * (1. - ρ_ygap) * (1. - ρ_smooth) # (1. - m[:ρ]) # This is the GDP part

    # Add MP shocks
    Γ0[eq[:eq_mp], endo[:rm_t]]     = -1.

    return Γ0, Γ1, C, Ψ, Π
end
