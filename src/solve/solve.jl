"""
```
solve(m::AbstractModel)
```

Driver to compute the model solution and augment transition matrices.

### Arguments
- `m`: the model object

### Return values
 - TTT, RRR, and CCC matrices of the state transition equation:
    ```
    S_t = TTT*S_{t-1} + RRR*ϵ_t + CCC
    ```
"""
function solve(m::AbstractModel)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = eqcond(m)

    # Solve model
    TTT_gensys, CCC_gensys, RRR_gensys, fmat, fwt, ywt, gev, eu, loose = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)
    if !((eu[1] == 1) & (eu[2] == 1))
        throw(GensysError("Gensys does not give existence"))
    end
    TTT_gensys = real(TTT_gensys)
    RRR_gensys = real(RRR_gensys)
    CCC_gensys = reshape(CCC_gensys, size(CCC_gensys, 1))

    # Augment states
    TTT, RRR, CCC = augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)

    return TTT, RRR, CCC
end



"""
```
GensysError <: Exception
```
A `GensysError` is thrown when Gensys does not give a unique solution, or no solution
exists. If a `GensysError`is thrown during Metropolis-Hastings, it is caught by `posterior`.
`posterior` then returns a value of `-Inf`, which Metropolis-Hastings always rejects.
### Fields
* `msg::ASCIIString`: Info message. Default = "Error in gensys." 
"""
type GensysError <: Exception
    msg::ASCIIString
end
GensysError() = GensysError("Error in gensys.")
Base.showerror(io::IO, ex::GensysError) = print(io, ex.msg)
