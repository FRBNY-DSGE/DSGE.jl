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

    Optim.optimize(fcn, x0,
                   method = Optim.NelderMead(parameters = parameters,
                                             initial_simplex = initial_simplex),
                   iterations = iterations, store_trace = store_trace, show_trace = show_trace,
                   extended_trace = extended_trace)
end

mutable struct MatlabSimplexer <: Optim.Simplexer
    a::Float64
    b::Float64
end
MatlabSimplexer(;a = 0.00025, b = 0.05) = MatlabSimplexer(a, b)

function Optim.simplexer(A::MatlabSimplexer, initial_x::Array{T, N}) where {T, N}
    n = length(initial_x)
    initial_simplex = Array{T, N}[initial_x for i = 1:n+1]
    for j = 1:n
        initial_simplex[j+1][j] += initial_simplex[j+1][j] == zero(T) ? A.b * initial_simplex[j+1][j] : A.a
    end
    initial_simplex
end
