function nelder_mead(fcn::Function,
                     x0::Array,
                     args...;
                     parameters            = Optim.AdaptiveParameters(),
                     initial_simplex       = Optim.AffineSimplexer(),
                     iterations::Int       = 1000,
                     store_trace::Bool     = false,
                     show_trace::Bool      = false,
                     extended_trace::Bool  = false,
                     kwargs...)

    # Modern Optim takes the method positionally and the run controls in an Options object
    # (the old `method = …` / loose `iterations = …` keyword API was removed).
    Optim.optimize(fcn, x0,
                   Optim.NelderMead(parameters = parameters, initial_simplex = initial_simplex),
                   Optim.Options(iterations = iterations, store_trace = store_trace,
                                 show_trace = show_trace, extended_trace = extended_trace))
end

mutable struct MatlabSimplexer <: Optim.Simplexer
    a::Float64
    b::Float64
end
MatlabSimplexer(;a = 0.00025, b = 0.05) = MatlabSimplexer(a, b)

function Optim.simplexer(A::MatlabSimplexer, initial_x::Array{T, N}) where {T, N}
    n = length(initial_x)
    # copy(initial_x) per vertex: `[initial_x for i=…]` aliases one array, so the loop
    # below would mutate every vertex at once and the simplex would be degenerate.
    initial_simplex = Array{T, N}[copy(initial_x) for i = 1:n+1]
    for j = 1:n
        # Match MATLAB fminsearch: a zero coordinate gets an absolute bump (A.a), a nonzero
        # one a relative bump (A.b * x_j). (The branches were previously swapped.)
        initial_simplex[j+1][j] += initial_simplex[j+1][j] == zero(T) ?
            A.a : A.b * initial_simplex[j+1][j]
    end
    initial_simplex
end
