"""
```
augment_states{T<:AbstractFloat}(m::AbstractDSGEModel, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T})
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
function augment_states(m::Model904, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where {T<:AbstractFloat}
    endo     = m.endogenous_states
    endo_new = m.endogenous_states_augmented
    exo      = m.exogenous_shocks

    n_endo = n_states(m)
    n_exo  = n_shocks_exogenous(m)
    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo)  == size(RRR)
    @assert n_endo           == length(CCC)

    # Initialize augmented matrices
    n_states_add = 7
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

    # Expected inflation
    TTT_aug[endo_new[:Et_π_t], 1:n_endo] = (TTT^2)[endo[:π_t], :]

    ### RRR modfications

    # Expected inflation
    RRR_aug[endo_new[:Et_π_t], :] = (TTT*RRR)[endo[:π_t], :]

    ### CCC Modifications

    # Expected inflation
    CCC_aug[endo_new[:Et_π_t]] = (CCC + TTT*CCC)[endo[:π_t]]

    return TTT_aug, RRR_aug, CCC_aug
end
