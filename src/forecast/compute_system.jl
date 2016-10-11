"""
`compute_system{T<:AbstractModel}(m::T)`

Solves the model `m` given the current parameter values for transition
matrices and computes measurement equation. Returns a `System` object.
"""
function compute_system{T<:AbstractModel}(m::T)

    # Solve model
    TTT, RRR, CCC = solve(m)
    transition_equation = Transition(TTT, RRR, CCC)
    
    # Solve measurement equation
    shocks = n_anticipated_shocks(m) > 0
    measurement_equation = measurement(m, TTT, RRR, CCC; shocks = shocks)
    
    return System(transition_equation, measurement_equation)
    
end
