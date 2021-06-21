"""
```
augment_states(m::AnSchorfheide, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T}) where {T<:AbstractFloat}
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
function augment_states(m::AnSchorfheide, TTT::Matrix{T}, RRR::Matrix{T},
                        CCC::Vector{T}; reg = 1) where {T<:AbstractFloat} # do not edit inputs
# this function needs to exist, but you probably don't need to edit it
# and can just leave it effectively an empty function (by immediately returning TTT, RRR, CCC)
    return TTT, RRR, CCC
end
