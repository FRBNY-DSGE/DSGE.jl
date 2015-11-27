using Compat, DSGE, HDF5

#=
doc"""
solve(m::AbstractModel)

### Parameters
- `m`: the model object

### Return values
 - TTT, RRR, and CCC matrices of the state transition equation:
              S_t = TTT*S_{t-1} + RRR*ϵ_t + CCC

## Description
Loads in the matrices that form the canonical representation of the equilibrium conditions by calling `eqcond`. Then calls `gensys` for the state-space representation of the model. Finally, calls `augment_states` to add growth rates to the observables. See documentation for each of these 3 functions for further details.
"""
=#
function solve(m::AbstractModel)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = eqcond(m) 
    
    # Solve model
    TTT_gensys, CCC_gensys, RRR_gensys = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
    TTT_gensys = real(TTT_gensys)
    RRR_gensys = real(RRR_gensys)
    CCC_gensys = reshape(CCC_gensys, length(CCC_gensys), 1)
    
    # Augment states
    TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)
    
    return TTT, RRR, CCC
end


#=
doc"""
augment_states{T<:AbstractFloat}(m::AbstractModel, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T})

### Parameters
-`m`: the m object
-`TTT`, `RRR`, and `CCC` are matrices of the state transition equation that is returned from `gensys`:
              S_t = TTT*S_{t-1} + RRR*ϵ_t + CCC

### Return values
- `TTT_aug`, `RRR_aug`, and `CCC_aug`, which extend the corresponding input matrices to include observables which are growth rates.

### Description
Some observables in the model are growth rates, which are calculated as a linear combination of a present and lagged state (which is not yet accounted for in the `TTT`, `RRR`,and `CCC` matrices). To improve the performance of `gensys`, these additional states are added after the model is solved. `augment_states` assigns an index to each lagged state, and extends the input `TTT`, `RRR`, and `CCC` matrices to accommodate the additional states and capture the lagged state value in the current state vector. `RRR` and `CCC` are mostly augmented with zeros.

The diagram below shows how `TTT` is extended to `TTT_aug`.

                TTT_aug
     (m.endogenous_states_postgensys
                  x
     m.endogenous_states_postgensys)                               
     ________________________ _________     
    |                        |         |
    |                        |         |
    |          TTT           |  lag    |
    | (m.endogenous_states   | states  |  
    |          x             |         |
    |  m.endogenous_states)  |         |
    |                        |         |
    |                        |         |
    |________________________|         |
    |                                  |
    |          lag                     |
    |         states                   |  
    |__________________________________|

"""                                
=#
function augment_states{T<:AbstractFloat}(m::Model990{T}, TTT::Matrix{T}, RRR::Matrix{T}, CCC::Matrix{T})
    endo = m.endogenous_states
    endo_addl = m.endogenous_states_postgensys
    exo = m.exogenous_shocks

    n_endo = n_states(m)
    n_exo = n_shocks_exogenous(m)
    @assert (n_endo, n_endo) == size(TTT)
    @assert (n_endo, n_exo) == size(RRR)
    @assert (n_endo, 1) == size(CCC)
    
    # Initialize augmented matrices
    numAdd = 12
    TTT_aug = zeros(n_endo + numAdd, n_endo + numAdd)
    TTT_aug[1:n_endo, 1:n_endo] = TTT
    RRR_aug = [RRR; zeros(numAdd, n_exo)]
    CCC_aug = [CCC; zeros(numAdd, 1)]

    ### TTT modifications

    # Track Lags
    TTT_aug[endo_addl[:y_t1], endo[:y_t]] = 1.0
    TTT_aug[endo_addl[:c_t1], endo[:c_t]] = 1.0
    TTT_aug[endo_addl[:i_t1], endo[:i_t]] = 1.0
    TTT_aug[endo_addl[:w_t1], endo[:w_t]] = 1.0
    TTT_aug[endo_addl[:π_t1], endo[:π_t]] = 1.0
    TTT_aug[endo_addl[:L_t1], endo[:L_t]]  = 1.0
    TTT_aug[endo_addl[:u_t1], endo[:u_t]] = 1.0

    # Expected inflation
    TTT_aug[endo_addl[:Et_π_t], 1:n_endo] = (TTT^2)[endo[:π_t], :]

    # The 8th column of the addition to TTT corresponds to "v_lr" which is set equal to
    # e_lr – measurement errors for the two real wage observables built in
    # as exogenous structural shocks.
    TTT_aug[endo_addl[:lr_t], endo_addl[:lr_t]] = m[:ρ_lr]
    TTT_aug[endo_addl[:tfp_t], endo_addl[:tfp_t]] = m[:ρ_tfp]
    TTT_aug[endo_addl[:e_gdpdef_t], endo_addl[:e_gdpdef_t]] = m[:ρ_gdpdef]
    TTT_aug[endo_addl[:e_pce_t], endo_addl[:e_pce_t]] = m[:ρ_pce]


    ### RRR modfications

    # Expected inflation
    RRR_aug[endo_addl[:Et_π_t], :] = (TTT*RRR)[endo[:π_t], :]

    # Measurement Error on long rate
    RRR_aug[endo_addl[:lr_t], exo[:lr_sh]] = 1.0

    # Measurement Error on TFP
    RRR_aug[endo_addl[:tfp_t], exo[:tfp_sh]] = 1.0

    # Measurement Error on GDP Deflator
    RRR_aug[endo_addl[:e_gdpdef_t], exo[:gdpdef_sh]] = 1.0

    # Measurement Error on Core PCE
    RRR_aug[endo_addl[:e_pce_t], exo[:pce_sh]] = 1.0



    ### CCC Modifications

    # Expected inflation
    CCC_aug[endo_addl[:Et_π_t], :] = (CCC + TTT*CCC)[endo[:π_t], :]

    return TTT_aug, RRR_aug, CCC_aug
end
