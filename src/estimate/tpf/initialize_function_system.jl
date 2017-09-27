"""
```
initialize_function_system{S<:AbstractFloat}(system::System{S})
```
### Inputs
- `system::System`: The output of compute_system(m), i.e. the matrix outputs from solving a given model, m.
### Output
- Returns the transition and measurement equations as functions,and the distributions of the shocks
and measurement error.
"""
function initialize_function_system{S<:AbstractFloat}(system::System{S})
    # Unpack system
    RRR    = system[:RRR]
    TTT    = system[:TTT]
    HH     = system[:EE] + system[:MM]*system[:QQ]*system[:MM]'
    DD     = system[:DD]
    ZZ     = system[:ZZ]
    QQ     = system[:QQ]
    EE     = system[:EE]
    MM     = system[:MM]
    sqrtS2 = RRR*get_chol(QQ)'

    @inline Φ(s_t1::Vector{S}, ε_t1::Vector{S}) = TTT*s_t1 + sqrtS2*ε_t1
    @inline Ψ(s_t1::Vector{S}, u_t1::Vector{S}) = ZZ*s_t1 + DD + u_t1

    F_ε = Distributions.MvNormal(zeros(size(QQ)[1]), eye(size(QQ)[1]))
    F_u = Distributions.MvNormal(zeros(size(HH)[1]), HH)

    return Φ, Ψ, F_ε, F_u
end
