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
function augment_states(m::RealBondMkup, TTT::Matrix{T}, TTT_jump::Matrix{T}, RRR::Matrix{T}, CCC::Vector{T}) where T<:AbstractFloat
    endo     = m.endogenous_states
    endo_unnorm = m.endogenous_states_unnormalized
    endo_new = m.model_states_augmented
    exo      = m.exogenous_shocks

    _n_model_states     = n_model_states(m)
    _n_states           = n_backward_looking_states(m)
    _n_jumps            = n_jumps(m)
    _n_observables      = n_observables(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    @assert (_n_model_states, _n_model_states) == size(TTT)
    @assert (_n_model_states, _n_shocks_exogenous)  == size(RRR)
    @assert (_n_model_states,)        == size(CCC)

    # Initialize augmented matrices
    n_states_add = length(endo_new)
    TTT_aug = zeros(_n_model_states + n_states_add, _n_model_states + n_states_add)
    TTT_aug[1:_n_model_states, 1:_n_model_states] = TTT
    RRR_aug = [RRR; zeros(n_states_add, _n_shocks_exogenous)]
    CCC_aug = [CCC; zeros(n_states_add)]

    #all of this is to compute GDPeqn or what on Marco's shee is Z^y
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
        construct_GDPfn_realbond(nx, ns, xgrid_total, sgrid,
                                 xwts, swts, γ, ν, abar, R, aborrow, μ, c, η, ell, χss)
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
    GDPeqn = GDPfn'*[eye(n_backward_looking_states_unnormalized(m)); gx2]*Qx'

    ### TTT modifications

    # Track Lags
            # Aug[:y_t1]                      States[:y_t]
   #= @show endo_new[:y_t1]
    @show size(TTT_aug)
    @show size(GDPeqn)=#
    TTT_aug[endo_new[:y_t1], 1:_n_states] = GDPeqn

    return TTT_aug, RRR_aug, CCC_aug
end
