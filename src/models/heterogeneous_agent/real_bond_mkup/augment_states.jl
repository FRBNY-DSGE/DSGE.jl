"""
```
augment_states{T<:AbstractFloat}(m::AbstractModel, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T})
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
function augment_states{T<:AbstractFloat}(m::RealBondMkup, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T})
    endo     = m.endogenous_states
    endo_unnorm = m.endogenous_states_unnormalized
    endo_new = m.model_states_augmented
    model_states_aug = m.model_states_augmented
    exo      = m.exogenous_shocks

    _n_model_states     = n_model_states(m)
    _n_states           = n_backward_looking_states(m)
    _n_jumps            = n_jumps(m)
    _n_observables      = n_observables(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    n_exo  = n_shocks_exogenous(m)
    @assert (_n_model_states, _n_model_states) == size(TTT)
    @assert (_n_model_states, n_exo)  == size(RRR)
    @assert (_n_model_states,)        == size(CCC)

    # Initialize augmented matrices
    n_states_add = _n_states #length(model_states_aug)
    TTT_aug = zeros(_n_model_states + n_states_add, _n_model_states + n_states_add)
    TTT_aug[1:_n_model_states, 1:_n_model_states] = TTT
    RRR_aug = [RRR; zeros(n_states_add, n_exo)]
    CCC_aug = [CCC; zeros(n_states_add)]

    ### TTT modifications

    # Track Lags
    @show endo
    @show endo_new
    @show _n_states
            # Aug[:y_t1]                      States[:y_t]
    TTT_aug[_n_model_states:size(TTT_aug, 1), 1:_n_states] = 1.0

    return TTT_aug, RRR_aug, CCC_aug
end
