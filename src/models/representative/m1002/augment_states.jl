"""
```
augment_states(m::AbstractDSGEModel, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T}) where {T<:AbstractFloat}
```

### Arguments

-`m`: the m object
-`TTT`, `RRR`, and `CCC`: matrices of the state transition equation

### Return values

- `TTT_aug`, `RRR_aug`, and `CCC_aug`: extend the corresponding input matrices to include
  jobservables which are growth rates.

### Description

Some observables in the model are growth rates, which are calculated as a linear combination
of a present and lagged state (which is not yet accounted for in the `TTT`, `RRR`,and `CCC`
matrices). To improve the performance of `gensys`, these additional states are added after
the model is solved. `augment_states` assigns an index to each lagged state, and extends the
input `TTT`, `RRR`, and `CCC` matrices to accommodate the additional states and capture the
lagged state value in the current state vector. `RRR` and `CCC` are mostly augmented with
zeros.

The diagram below shows how `TTT` is extended to `TTT_aug`.

                TTT_aug
     (m.endogenous_states_additional
                  x
     m.endogenous_states_additional)
     _________________________________
    |                     |           |
    |          TTT        | endog_    |
    | (endogenous_states  | states_   |
    |          x          | augmented |
    |  endogenous_states) |           |
    |_____________________|           |
    |                                 |
    |    endogenous_states_augmented  |
    |_________________________________|

"""
function augment_states(m::Model1002, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T};
                        regime_switching::Bool = false, reg::Int = 1) where {T<:AbstractFloat}
    endo     = m.endogenous_states
    endo_new = m.endogenous_states_augmented
    exo      = m.exogenous_shocks

    n_endo = n_states(m)
    n_exo  = n_shocks_exogenous(m)
    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo)  == size(RRR)
    @assert (n_endo,)        == size(CCC)

    # Initialize augmented matrices
    n_states_add = length(endo_new)

    TTT_aug = zeros(n_endo + n_states_add, n_endo + n_states_add)
    TTT_aug[1:n_endo, 1:n_endo] = TTT
    RRR_aug = [RRR; zeros(n_states_add, n_exo)]
    CCC_aug = [CCC; zeros(n_states_add)]

    for para in m.parameters
        if !isempty(para.regimes)
            if (haskey(get_settings(m), :model2para_regime) ? haskey(get_setting(m, :model2para_regime), para.key) : false)
                ModelConstructors.toggle_regime!(para, reg, get_setting(m, :model2para_regime)[para.key])
            else
                ModelConstructors.toggle_regime!(para, reg)
            end
        end
    end

    ### TTT modifications

    # Track Lags
    TTT_aug[endo_new[:y_t1],     endo[:y_t]] = 1.0
    TTT_aug[endo_new[:c_t1],     endo[:c_t]] = 1.0
    TTT_aug[endo_new[:i_t1],     endo[:i_t]] = 1.0
    TTT_aug[endo_new[:w_t1],     endo[:w_t]] = 1.0
    TTT_aug[endo_new[:π_t1_dup], endo[:π_t]] = 1.0
    TTT_aug[endo_new[:L_t1],     endo[:L_t]] = 1.0
    TTT_aug[endo_new[:u_t1],     endo[:u_t]] = 1.0
    TTT_aug[endo_new[:e_gdp_t1], endo_new[:e_gdp_t]] = 1.0
    TTT_aug[endo_new[:e_gdi_t1], endo_new[:e_gdi_t]] = 1.0
    if subspec(m) in ["ss14", "ss15", "ss16", "ss18", "ss19"]
        TTT_aug[endo_new[:e_tfp_t1], endo_new[:e_tfp_t]] = 1.0
    end

    if get_setting(m, :add_laborproductivity_measurement)
        TTT_aug[endo_new[:cum_z_t],  endo[:z_t]]         = 1.0
        TTT_aug[endo_new[:cum_z_t],  endo_new[:cum_z_t]] = 1.0
        CCC_aug[endo_new[:cum_z_t]]                      = 100. * (exp(m[:z_star]) - 1.)
    end

    if get_setting(m, :add_nominalgdp_level)
        TTT_aug[endo_new[:cum_y_t],  endo[:y_t]]         = 1.0
        TTT_aug[endo_new[:cum_y_t],  endo_new[:y_t1]]    = -1.0
        TTT_aug[endo_new[:cum_y_t],  endo_new[:cum_y_t]] = 1.0

        TTT_aug[endo_new[:cum_z_t],  endo[:z_t]]         = 1.0
        TTT_aug[endo_new[:cum_z_t],  endo_new[:cum_z_t]] = 1.0
        CCC_aug[endo_new[:cum_z_t]]                      = 100. * (exp(m[:z_star]) - 1.)

        TTT_aug[endo_new[:cum_π_t],  endo[:π_t]]         = 1.0
        TTT_aug[endo_new[:cum_π_t],  endo_new[:cum_π_t]] = 1.0
        CCC_aug[endo_new[:cum_π_t]]                      = 100. * (m[:π_star] - 1.)

        TTT_aug[endo_new[:cum_e_gdp_t],  endo_new[:e_gdp_t]]     = 1.0
        TTT_aug[endo_new[:cum_e_gdp_t],  endo_new[:e_gdp_t1]]    = -m[:me_level]
        TTT_aug[endo_new[:cum_e_gdp_t],  endo_new[:cum_e_gdp_t]] = 1.0

        # Consumption
        TTT_aug[endo_new[:cum_c_t],  endo[:c_t]]         = 1.0
        TTT_aug[endo_new[:cum_c_t],  endo_new[:c_t1]]    = -1.0
        TTT_aug[endo_new[:cum_c_t],  endo_new[:cum_c_t]] = 1.0

        # Flexible Consumption
        TTT_aug[endo_new[:cum_c_f_t],  endo[:c_f_t]]         = 1.0
        TTT_aug[endo_new[:cum_c_f_t],  endo_new[:c_f_t1]]    = -1.0
        TTT_aug[endo_new[:cum_c_f_t],  endo_new[:cum_c_f_t]] = 1.0
        TTT_aug[endo_new[:c_f_t1],     endo[:c_f_t]]         = 1.0

        # Investment
        TTT_aug[endo_new[:cum_i_t],  endo[:i_t]]         = 1.0
        TTT_aug[endo_new[:cum_i_t],  endo_new[:i_t1]]    = -1.0
        TTT_aug[endo_new[:cum_i_t],  endo_new[:cum_i_t]] = 1.0

        # Flexible Investment
        TTT_aug[endo_new[:cum_i_f_t],  endo[:i_f_t]]         = 1.0
        TTT_aug[endo_new[:cum_i_f_t],  endo_new[:i_f_t1]]    = -1.0
        TTT_aug[endo_new[:cum_i_f_t],  endo_new[:cum_i_f_t]] = 1.0
        TTT_aug[endo_new[:i_f_t1],     endo[:i_f_t]]         = 1.0
    end

    if get_setting(m, :add_flexible_price_growth)
        # Track additional lags
        TTT_aug[endo_new[:y_f_t1],     endo[:y_f_t]]         = 1.0
        TTT_aug[endo_new[:c_f_t1],     endo[:c_f_t]]         = 1.0
        TTT_aug[endo_new[:i_f_t1],     endo[:i_f_t]]         = 1.0
    end

    if get_setting(m, :add_cumulative)
        # Output Gap
        TTT_aug[endo_new[:cum_y_t],  endo[:y_t]]         = 1.0
        TTT_aug[endo_new[:cum_y_t],  endo_new[:y_t1]]    = -1.0
        TTT_aug[endo_new[:cum_y_t],  endo_new[:cum_y_t]] = 1.0

        TTT_aug[endo_new[:cum_y_f_t],  endo[:y_f_t]]         = 1.0
        TTT_aug[endo_new[:cum_y_f_t],  endo_new[:y_f_t1]]    = -1.0
        TTT_aug[endo_new[:cum_y_f_t],  endo_new[:cum_y_f_t]] = 1.0
        TTT_aug[endo_new[:y_f_t1],     endo[:y_f_t]]         = 1.0

        # Remaining accumulations for GDP, Flexible GDP, Technology
        TTT_aug[endo_new[:cum_z_t],  endo[:z_t]]         = 1.0
        TTT_aug[endo_new[:cum_z_t],  endo_new[:cum_z_t]] = 1.0
        CCC_aug[endo_new[:cum_z_t]]                      = 100. * (exp(m[:z_star]) - 1.)

        TTT_aug[endo_new[:cum_e_gdp_t],  endo_new[:e_gdp_t]]     = 1.0
        TTT_aug[endo_new[:cum_e_gdp_t],  endo_new[:e_gdp_t1]]    = -m[:me_level]
        TTT_aug[endo_new[:cum_e_gdp_t],  endo_new[:cum_e_gdp_t]] = 1.0

        # Consumption
        TTT_aug[endo_new[:cum_c_t],  endo[:c_t]]         = 1.0
        TTT_aug[endo_new[:cum_c_t],  endo_new[:c_t1]]    = -1.0
        TTT_aug[endo_new[:cum_c_t],  endo_new[:cum_c_t]] = 1.0

        # Flexible Consumption
        TTT_aug[endo_new[:cum_c_f_t],  endo[:c_f_t]]         = 1.0
        TTT_aug[endo_new[:cum_c_f_t],  endo_new[:c_f_t1]]    = -1.0
        TTT_aug[endo_new[:cum_c_f_t],  endo_new[:cum_c_f_t]] = 1.0
        TTT_aug[endo_new[:c_f_t1],     endo[:c_f_t]]         = 1.0

        # Investment
        TTT_aug[endo_new[:cum_i_t],  endo[:i_t]]         = 1.0
        TTT_aug[endo_new[:cum_i_t],  endo_new[:i_t1]]    = -1.0
        TTT_aug[endo_new[:cum_i_t],  endo_new[:cum_i_t]] = 1.0

        # Flexible Investment
        TTT_aug[endo_new[:cum_i_f_t],  endo[:i_f_t]]         = 1.0
        TTT_aug[endo_new[:cum_i_f_t],  endo_new[:i_f_t1]]    = -1.0
        TTT_aug[endo_new[:cum_i_f_t],  endo_new[:cum_i_f_t]] = 1.0
        TTT_aug[endo_new[:i_f_t1],     endo[:i_f_t]]         = 1.0
    end

    if get_setting(m, :add_flexible_price_growth)
        # Track additional lags
        TTT_aug[endo_new[:y_f_t1],     endo[:y_f_t]]         = 1.0
        TTT_aug[endo_new[:c_f_t1],     endo[:c_f_t]]         = 1.0
        TTT_aug[endo_new[:i_f_t1],     endo[:i_f_t]]         = 1.0
    end

    # Expected inflation
    # TTT_aug[endo_new[:Et_π_t], 1:n_endo] = (TTT^2)[endo[:π_t], :]
    TTT_aug[endo_new[:Et_π_t], 1:n_endo] = view(TTT^2, endo[:π_t], :)

    # The 8th column of the addition to TTT corresponds to "v_lr" which is set equal to
    # e_lr – measurement errors for the two real wage observables built in
    # as exogenous structural shocks.
    TTT_aug[endo_new[:e_lr_t], endo_new[:e_lr_t]]           = m[:ρ_lr]
    TTT_aug[endo_new[:e_tfp_t], endo_new[:e_tfp_t]]         = m[:ρ_tfp]
    TTT_aug[endo_new[:e_gdpdef_t], endo_new[:e_gdpdef_t]]   = m[:ρ_gdpdef]
    TTT_aug[endo_new[:e_corepce_t], endo_new[:e_corepce_t]] = m[:ρ_corepce]
    TTT_aug[endo_new[:e_gdp_t], endo_new[:e_gdp_t]]         = m[:ρ_gdp]
    TTT_aug[endo_new[:e_gdi_t], endo_new[:e_gdi_t]]         = m[:ρ_gdi]

    if subspec(m) in ["ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78"]
        # COVID counterparts to measurement errors
        TTT_aug[endo_new[:e_lr_t], endo_new[:e_lr_covid_t]]   = m[:ρ_lr_covid]
        TTT_aug[endo_new[:e_tfp_t], endo_new[:e_tfp_covid_t]] = m[:ρ_tfp_covid]
        TTT_aug[endo_new[:e_gdp_t], endo_new[:e_gdp_covid_t]] = m[:ρ_gdp_covid]
        TTT_aug[endo_new[:e_gdi_t], endo_new[:e_gdi_covid_t]] = m[:ρ_gdi_covid]
    end

    if haskey(get_settings(m), :add_iid_cond_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_cond_obs_gdp_meas_err) : false
        TTT_aug[endo_new[:e_condgdp_t], endo_new[:e_condgdp_t]] = m[:ρ_condgdp]
    end

    if haskey(get_settings(m), :add_iid_anticipated_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_anticipated_obs_gdp_meas_err) : false
        TTT_aug[endo_new[:e_gdpexp_t], endo_new[:e_gdpexp_t]] = m[:ρ_gdpexp]
    end

    if haskey(get_settings(m), :add_iid_cond_obs_corepce_meas_err) ?
        get_setting(m, :add_iid_cond_obs_corepce_meas_err) : false
        TTT_aug[endo_new[:e_condcorepce_t], endo_new[:e_condcorepce_t]] = m[:ρ_condcorepce]
    end

    # Fundamental inflation
    if subspec(m) in ["ss13", "ss14", "ss15", "ss16", "ss17", "ss18", "ss19"]
        betabar = exp((1-m[:σ_c] ) * m[:z_star]) * m[:β]
        κ = ((1 - m[:ζ_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))*
             (1 - m[:ζ_p]))/(m[:ζ_p]*((m[:Φ]- 1)*m[:ϵ_p] + 1))/
        (1 + m[:ι_p]*m[:β]*exp((1 - m[:σ_c])*m[:z_star]))
        κcoef = κ * (1 + m[:ι_p] * betabar)
        sₜ∞ = inv(Matrix{Float64}(I,size(TTT)...) - betabar * TTT) # infinite horizon expectation of state vector

        TTT_aug[endo_new[:Sinf_t], 1:n_endo] = sₜ∞[endo[:mc_t],:]
        TTT_aug[endo_new[:πtil_t], endo_new[:Sinf_t]] = κcoef
        TTT_aug[endo_new[:πtil_t], endo_new[:πtil_t1]] = m[:ι_p]
        TTT_aug[endo_new[:πtil_t1], endo_new[:πtil_t]] = 1.
    end

    ### RRR modfications

    # Expected inflation
    # RRR_aug[endo_new[:Et_π_t], :] = (TTT*RRR)[endo[:π_t], :]
    RRR_aug[endo_new[:Et_π_t], :] = view(TTT * RRR, endo[:π_t], :)

    # Measurement Error on long rate
    RRR_aug[endo_new[:e_lr_t], exo[:lr_sh]] = 1.0

    # Measurement Error on TFP
    RRR_aug[endo_new[:e_tfp_t], exo[:tfp_sh]] = 1.0

    # Measurement Error on GDP Deflator
    RRR_aug[endo_new[:e_gdpdef_t], exo[:gdpdef_sh]] = 1.0

    # Measurement Error on Core PCE
    RRR_aug[endo_new[:e_corepce_t], exo[:corepce_sh]] = 1.0

    # Measurement Error on GDP and GDI
    RRR_aug[endo_new[:e_gdp_t], exo[:gdp_sh]] = 1.0
    RRR_aug[endo_new[:e_gdp_t], exo[:gdi_sh]] = m[:ρ_gdpvar] * m[:σ_gdp] ^ 2

    RRR_aug[endo_new[:e_gdi_t], exo[:gdi_sh]] = 1.0

    if subspec(m) in ["ss67", "ss68", "ss69", "ss70", "ss71", "ss72", "ss73", "ss74", "ss75", "ss76", "ss77", "ss78"]
        # COVID counterparts to measurement errors
        RRR_aug[endo_new[:e_lr_t], exo[:lr_covid_sh]]   = 1.0
        RRR_aug[endo_new[:e_tfp_t], exo[:tfp_covid_sh]] = 1.0

        RRR_aug[endo_new[:e_gdp_t], exo[:gdp_covid_sh]] = 1.0
        RRR_aug[endo_new[:e_gdp_t], exo[:gdi_covid_sh]] = m[:ρ_gdpvar_covid] * m[:σ_gdp_covid] ^ 2

        RRR_aug[endo_new[:e_gdi_t], exo[:gdi_covid_sh]] = 1.0
    end

    if haskey(get_settings(m), :add_iid_cond_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_cond_obs_gdp_meas_err) : false
        RRR_aug[endo_new[:e_condgdp_t], exo[:condgdp_sh]] = 1.0
    end

    if haskey(get_settings(m), :add_iid_anticipated_obs_gdp_meas_err) ?
        get_setting(m, :add_iid_anticipated_obs_gdp_meas_err) : false
        RRR_aug[endo_new[:e_gdpexp_t], exo[:gdpexp_sh]] = 1.0
    end

    if haskey(get_settings(m), :add_iid_cond_obs_corepce_meas_err) ?
        get_setting(m, :add_iid_cond_obs_corepce_meas_err) : false
        RRR_aug[endo_new[:e_condcorepce_t], exo[:condcorepce_sh]] = 1.0
    end

    ### CCC Modifications

    # Expected inflation
    CCC_aug[endo_new[:Et_π_t]] = (CCC + TTT*CCC)[endo[:π_t]] # note that currently this term is not used anywhere


    ### Ensure parameter values are in regime 1
    for para in m.parameters
        if !isempty(para.regimes)
            ModelConstructors.toggle_regime!(para, 1)
        end
    end

    return TTT_aug, RRR_aug, CCC_aug
end
