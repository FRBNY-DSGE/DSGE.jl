# creates the new type ReductionData to hold information on the type of reduction(s) being performed
# Fields
# reduce_distribution: True if you want to reduce distribution via Krylov reduction
# reduce_v: True if you want to use spline basis reduction on the value function
# krylov_dim: Desired dimension of Krylov subspace. Initialized to 0 if unused.
# knots_dict: Dictionaries of knot points for each dimension along which value function is approximated.
#             Eempty dict if unused.
# ss_array: Array giving the state space along which the spline basis approximation is done.
#           Empty array if unused.
# n_prior: Number of grid points stacked prior to reduction dimension. 0 if unused.
# n_post: Number of grid points stacked after reduction dimension. 0 if unused.

 struct ReductionData
     reduce_distribution::Bool
     reduce_v::Bool
     krylov_dim::Int64
     knots_dict::Dict{Int64, Vector{Float64}}
     ss_array::Array{Float64}
     n_prior::Int64
     n_post::Int64
end

function ReductionData(reduce_distribution::Bool, reduce_v::Bool)
    if !reduce_distribution && !reduce_v
        ReductionData(reduce_distribution, reduce_v, 0 , Dict{Int64, Vector{Float64}}(), Array{Float64}(), 0, 0)
    else
        error("If not both false, must specify additional fields to perform reductions.")
    end
end

function ReductionData(reduce_distribution::Bool, reduce_v::Bool, krylov_dim::Int64)
    if reduce_distribution && !reduce_v
        ReductionData(reduce_distribution, reduce_v, krylov_dim , Dict{Int64, Vector{Float64}}(), Array{Float64}(), 0, 0)
    else
        error("Input Error: this method allows for just distribution reduction.")
    end
end

function ReductionData(reduce_distribution::Bool, reduce_v::Bool, knots_dict::Dict{Int64, Vector{Float64}},
                       ss_array::Array{Float64}, n_prior::Int64, n_post::Int64)
    if reduce_v && !reduce_distribution
        ReductionData(reduce_distribution, reduce_v, 0 , knots_dict, ss_array, n_prior, n_post)
    else
        error("Input Error: this method allows for just value function reduction.")
    end
end