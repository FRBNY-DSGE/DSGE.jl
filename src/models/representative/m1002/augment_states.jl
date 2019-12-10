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
function augment_states(m::Model1002, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
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

    # Expected inflation
    TTT_aug[endo_new[:Et_π_t], 1:n_endo] = (TTT^2)[endo[:π_t], :]

    # The 8th column of the addition to TTT corresponds to "v_lr" which is set equal to
    # e_lr – measurement errors for the two real wage observables built in
    # as exogenous structural shocks.
    TTT_aug[endo_new[:e_lr_t], endo_new[:e_lr_t]]               = m[:ρ_lr]
    TTT_aug[endo_new[:e_tfp_t], endo_new[:e_tfp_t]]             = m[:ρ_tfp]
    TTT_aug[endo_new[:e_gdpdef_t], endo_new[:e_gdpdef_t]]   = m[:ρ_gdpdef]
    TTT_aug[endo_new[:e_corepce_t], endo_new[:e_corepce_t]] = m[:ρ_corepce]
    TTT_aug[endo_new[:e_gdp_t], endo_new[:e_gdp_t]]         = m[:ρ_gdp]
    TTT_aug[endo_new[:e_gdi_t], endo_new[:e_gdi_t]]         = m[:ρ_gdi]

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
    RRR_aug[endo_new[:Et_π_t], :] = (TTT*RRR)[endo[:π_t], :]

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

    ### CCC Modifications

    # Expected inflation
    CCC_aug[endo_new[:Et_π_t]] = (CCC + TTT*CCC)[endo[:π_t]]

    return TTT_aug, RRR_aug, CCC_aug
end
