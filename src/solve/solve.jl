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


