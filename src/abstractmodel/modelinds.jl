# A type that bundles together the model indices functions from models/m$(spec)/modelinds.jl, as well as the number of states, etc.
type ModelInds
    endo::Function
    exo::Function
    ex::Function
    eq::Function

    # These are used to give the dimensions of the matrices in eqconds.jl
    n_states::Int64 # number of endogenous states
    n_exo::Int64 # number of exogenous (iid) shocks
    n_ex::Int64 # number of rational expectations errors
    n_eq::Int64
end



# Takes as input the four functions defined in models/m$(spec)/modelinds.jl and returns a ModelInds object
function ModelInds(make_endo::Function, make_exo::Function, make_ex::Function, make_eq::Function)
    endo, n_states = make_endo()
    exo, n_exo = make_exo()
    ex, n_ex = make_ex() # so named to avoid a namespace collision with Base.exp()
    eq, n_eq = make_eq()

    return ModelInds(endo, exo, ex, eq, n_states, n_exo, n_ex, n_eq)
end
