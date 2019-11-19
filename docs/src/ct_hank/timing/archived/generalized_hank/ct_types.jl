```
This file holds several new data types to help users solve continuous time models,
particularly HANK style models. This file defines two types:
Grid, ValueFunction
```

# The following code creates a new data type Grid which stores information
# and associated methods related to state space grids
#
# Fields
# state_space: A dictionary whose keys are dimensions and whose values
#              are vectors holding the state space
# state_space_length: vector holding size of each dimension
# ss_array: keys -> dimension, values -> array of a dimension repeated over other dimensions
# df: keys -> dimension, values -> forward difference
# db: keys -> dimension, values -> backward difference
# df_grid: keys -> dimension, values -> repeated forward difference matrix
# db_grid: keys -> dimension, values -> repeated backward difference matrix
# delta: keys -> dimension, values -> lengths of overlapping intervals to allow integration
# deltamat: keys -> dimension, values -> repetition of delta across other dimensions
# auxiliary: keys -> dimension, values -> any auxiliary information for a given dimension
struct Grid
    state_space::Dict{Int64, Vector{Float64}}
    state_space_length::Vector{Int64}
    ss_array::Dict{Int64, Array{Float64}}
    df::Dict{Int64, Vector{Float64}}
    db::Dict{Int64, Vector{Float64}}
    df_grid::Dict{Int64, Array{Float64}}
    db_grid::Dict{Int64, Array{Float64}}
    integration_vec::Dict{Int64, Vector{Float64}}
    integration_grid::Dict{Int64, Array{Float64}}
    auxiliary::Dict{Int64, Dict{Symbol, Any}}
end

# Automatically fills in all the fields or does not fill in any, leaving the values
# for the user to define.
function Grid(state_space::Dict{Int64, Vector{Float64}}; fill::Bool = true)
    length_vec = map(i -> length(state_space[i]), 1:length(keys(state_space)))
    grid = Grid(state_space, length_vec, Dict{Int64, Array{Float64}}(), Dict{Int64, Vector{Float64}}(),
         Dict{Int64, Vector{Float64}}(), Dict{Int64, Array{Float64}}(), Dict{Int64, Array{Float64}}(),
         Dict{Int64, Vector{Float64}}(), Dict{Int64, Array{Float64}}(),
         Dict{Int64, Dict{Symbol, Any}}())
    if fill
        set_ss_array(grid)
        set_df(grid)
        set_db(grid)
        set_df_grid(grid; from_scratch = false)
        set_db_grid(grid; from_scratch = false)
        set_int_vec(grid; from_scratch = false)
        set_int_grid(grid; from_scratch = false)
    end
    return grid
end

# Computes ss_array based on state_space
function set_ss_array(grid::Grid; dims::Vector{Int64} = collect(1:length(keys(grid.state_space))))
    if length(dims) < 1
        error("No dimensions of state space entered.")
    end
    for dim in dims
        grid.ss_array[dim] = repeat(grid.state_space[dim],
                                        outer = hcat(1, grid.state_space_length[1:end .!= dim]'))
        if dim != 1
            permute_vector = [collect(2:dim)' 1 collect(dim + 1:length(dims))']
            grid.ss_array[dim] = permutedims(grid.ss_array[dim], permute_vector)
        end
    end
end

# Computes df
function set_df(grid::Grid; dims::Vector{Int64} = collect(1:length(keys(grid.state_space))))
    if length(dims) < 1
        error("No dimensions of state space entered.")
    end
    for dim in dims
        grid.df[dim] = similar(grid.state_space[dim])
        grid.df[dim][1:end - 1] = grid.state_space[dim][2:end] - grid.state_space[dim][1:end - 1]
        grid.df[dim][end] = grid.df[dim][end - 1]
    end
end

# Computes db
function set_db(grid::Grid; dims::Vector{Int64} = collect(1:length(keys(grid.state_space))))
    if length(dims) < 1
        error("No dimensions of state space entered.")
    end
    for dim in dims
        grid.db[dim] = similar(grid.state_space[dim])
        grid.db[dim][2:end] = grid.state_space[dim][2:end] - grid.state_space[dim][1:end - 1]
        grid.db[dim][1] = grid.db[dim][2]
    end
end

# Computes df_grid, if from_scratch is true, then we assume df was not computed yet
function set_df_grid(grid::Grid; from_scratch::Bool = true, save::Bool = false, dims::Vector{Int64} = collect(1:length(keys(grid.state_space))))
    if length(dims) < 1
        error("No dimensions of state space entered.")
    end
    if from_scratch
        if save
            df = grid.df
        else
            df = Dict{Int64, Vector{Float64}}()
        end
        for dim in dims
            df[dim] = similar(grid.state_space[dim])
            df[dim][1:end - 1] = grid.state_space[dim][2:end] - grid.state_space[dim][1:end - 1]
            df[dim][end] = df[dim][end - 1]
        end
    else
        df = grid.df
    end

    # Create df_grid matrices
    for dim in dims
        grid.df_grid[dim] = repeat(df[dim], outer = hcat(1, grid.state_space_length[1:end .!= dim]'))
        if dim != 1
            permute_vector = [collect(2:dim)' 1 collect(dim + 1:length(dims))']
            grid.df_grid[dim] = permutedims(grid.df_grid[dim], permute_vector)
        end
    end
end

# Computes db_grid
function set_db_grid(grid::Grid; from_scratch::Bool = false, save::Bool = false, dims::Vector{Int64} = collect(1:length(keys(grid.state_space))))
    if length(dims) < 1
        error("No dimensions of state space entered.")
    end
    if from_scratch
        if save
            db = grid.db
        else
            db = Dict{Int64, Vector{Float64}}()
        end
        for dim in dims
            db[dim] = similar(grid.state_space[dim])
            db[dim][2:end] = grid.state_space[dim][2:end] - grid.state_space[dim][1:end - 1]
            db[dim][1] = db[dim][2]
        end
    else
        db = grid.db
    end

    # Create db_grid matrices
    for dim in dims
        grid.db_grid[dim] = repeat(db[dim], outer = hcat(1, grid.state_space_length[1:end .!= dim]'))
        if dim != 1
            permute_vector = [collect(2:dim)' 1 collect(dim + 1:length(dims))']
            grid.db_grid[dim] = permutedims(grid.db_grid[dim], permute_vector)
        end
    end
end

# Computes integration_vec
function set_int_vec(grid::Grid; from_scratch::Bool = true, save::Bool = false, dims::Vector{Int64} = collect(1:length(keys(grid.state_space))))
    if length(dims) < 1
        error("No dimensions of state space entered.")
    end
    if from_scratch
        if save
            df = grid.df
        else
            df = Dict{Int64, Vector{Float64}}()
        end
        for dim in dims
            df[dim] = similar(grid.state_space[dim])
            df[dim][1:end - 1] = grid.state_space[dim][2:end] - grid.state_space[dim][1:end - 1]
            df[dim][end] = df[dim][end - 1]
        end
    else
        df = grid.df
    end

    # Create vectors
    for dim in dims
        grid.integration_vec[dim] = similar(grid.state_space[dim])
        grid.integration_vec[dim][1] = 0.5 * df[dim][1]
        grid.integration_vec[dim][2:end - 1] = 0.5 * (df[dim][1:end - 2] + df[dim][2:end - 1])
        grid.integration_vec[dim][end] = 0.5 * df[dim][end - 1]
    end
end

# Computes integration_grid
function set_int_grid(grid::Grid; from_scratch::Bool = true, save::Bool = false, dims::Vector{Int64} = collect(1:length(keys(grid.state_space))))
    if length(dims) < 1
        error("No dimensions of state space entered.")
    end
    if from_scratch
        if save
            int_vec = grid.integration_vec
p        else
            int_vec = Dict{Int64, Vector{Float64}}
        end
        df = Dict{Int64, Vector{Float64}}
        for dim in dims
            df[dim] = similar(grid.state_space[dim])
            df[dim][1:end - 1] = grid.state_space[dim][2:end] - grid.state_space[dim][1:end - 1]
            df[dim][end] = df[dim][end - 1]
        end
        for dim in dims
            int_vec[dim] = similar(grid.state_space[dim])
            int_vec[dim][1] = 0.5 * df[dim][1]
            int_vec[dim][2:end - 1] = 0.5 * (df[dim][1:end - 2] + df[dim][2:end - 1])
            int_vec[dim][end] = 0.5 * df[dim][end - 1]
        end
    else
        int_vec = grid.integration_vec
    end

    # Compute grid
    for dim in dims
        grid.integration_grid[dim] = repeat(int_vec[dim], outer = hcat(1, grid.state_space_length[1:end .!= dim]'))
        if dim != 1
            permute_vector = [collect(2:dim)' 1 collect(dim + 1:length(dims))']
            grid.integration_grid[dim] = permutedims(grid.integration_grid[dim], permute_vector)
        end
    end
end

# The following code creates ValueFunction,
# which holds data about the value function
# over the state space defined by a Grid type
# object and automatically computes
# relevant objects like forward difference matrices.
#
# Fields
# grid: Holds a Grid object defining the state space
# V: value function over state space
# Vf: forward difference matrices in each dimension's direction
# Vb: backward difference matrices in each dimension's direction

# struct ValueFunction
#     grid::Grid
#     V::Array{Float64}
#     Vf::Dict{Int64, Array{Float64}}
#     Vb::Dict{Int64, Array{Float64}}
# end

# # Main constructor
# function ValueFunction(grid::Grid, V::Array{Float64}; fill::Bool = true)
#     vf = ValueFunction(grid, V, Dict{Int64, Array{Float64}}(), Dict{Int64, Array{Float64}}())
#     if fill
#         set_Vf(vf)
#         set_Vb(vf)
#     end
# end

# # Computes forward differences except boundary conditions
# function set_Vf(val_fnct::ValueFunction; dims::Vector{Int64} = collect(1:length(keys(val_fnct.grid.state_space))))
#     for dim in dims
#         # Initialize forward difference matrix
#         vf = zeros(val_fnct.V)
#         if dim > 1
#             V_copy = copy(val_fnct.V)

#             # permute dimensions
#             permute_vec = dims
#             permute_vec[1] = dim
#             permute_vec[dim] = 1
#             V_copy = permutedims(V_copy, permute_vec)
#         end
#         cds = size(V_copy, 1) # cds for current dimension size
#         for i = 1:prod(grid.state_space_length[1:end .!= dim])
#             vf[1 + (i-1)*cds:cds - 1 + (i-1)*cds] = (V_copy[2 + (i-1)*cds:cds + (i-1)*cds] -
#                                                      V_copy[1 + (i-1)*cds:cds - 1 + (i-1)*cds]) ./
#                                                      val_fnct.grid.df[dim][1 + (i-1)*cds:cds - 1 + (i-1)*cds]
#         end
#         if dim > 1
#             V.vf[dim] = permutedims(vf, permute_vec)
#             permute_vec[1] = 1; permute_vec[dim] = dim
#         else
#             V.vf[dim] = vf
#         end
#     end
# end

# # Computes backward differences except boundary conditions
# function set_Vb(val_fnct::ValueFunction; dims::Vector{Int64} = collect(1:length(keys(vf.grid.state_space))))
#     for dim in dims
#         # Initialize forward difference matrix
#         vb = zeros(val_fnct.V)
#         if dim > 1
#             V_copy = copy(val_fnct.V)

#             # permute dimensions
#             permute_vec = dims
#             permute_vec[1] = dim
#             permute_vec[dim] = 1
#             V_copy = permutedims(V_copy, permute_vec)
#         end
#         cds = size(V_copy, 1) # cds for current dimension size
#         for i = 1:prod(grid.state_space_length[1:end .!= dim])
#             vb[2 + (i-1)*cds:cds + (i-1)*cds] = (V_copy[2 + (i-1)*cds:cds + (i-1)*cds] -
#                                                      V_copy[1 + (i-1)*cds:cds - 1 + (i-1)*cds]) ./
#                                                      val_fnct.grid.db[dim][2 + (i-1)*cds:cds + (i-1)*cds]
#         end
#         if dim > 1
#             V.vb[dim]= permutedims(vb, permute_vec)
#             permute_vec[1] = 1
#             permute_vec[dim] = dim
#         else
#             V.vb[dim] = vb
#         end
#     end
# end

# Allows user to specify forward difference boundary conditions along one dimension.
# Assumes either a stacked vector for the boundary conditions
# or an array that has the given dimension removed.
# function set_Vf_bndy(val_fnct::ValueFunction, bndy::Array{Float64}, dim::Int64)
#     # check that bndy is either stacked or has the right dimensions
#     prod_dims = prod(val_fnct.state_space_length[1:end .!= dim])
#     if !( (typeof(bndy) == Vector{Float64} && length(bndy) == prod_dims) ||
#           size(bndy) == tuple(val_fnct.state_space_length[1:end .!= dim]...))
#         error("Provided boundary conditions do not have the right dimensions.")
#     end

#     if dim > 1
#         Vf_copy = copy(val_fnct.Vf)

#         # permute dimensions
#         permute_vec = copy(val_fnct.state_space_length)
#         permute_vec[1] = dim
#         permute_vec[dim] = 1
#         Vf_copy = permutedims(Vf_copy, permute_vec)
#     else
#         Vf_copy = val_fnct.Vf
#     end
#     cds = size(Vf_copy, 1) # cds for current dimension size
#     for i = 1:prod(grid.state_space_length[1:end .!= dim])
#         Vf_copy[cds + (i-1)*cds] = bndy[i]
#     end
#     if dim > 1
#         val_fnct.Vf = permutedims(Vf_copy, permute_vec)
#     end
# end

# # Allows user to specify backward difference boundary conditions
# function set_Vf_bndy(val_fnct::ValueFunction, bndy::Array{Float64}, dim::Int64)
#     # check that bndy is either stacked or has the right dimensions
#     prod_dims = prod(val_fnct.state_space_length[1:end .!= dim])
#     if !( (typeof(bndy) == Vector{Float64} && length(bndy) == prod_dims) ||
#           size(bndy) == tuple(val_fnct.state_space_length[1:end .!= dim]...))
#         error("Provided boundary conditions do not have the right dimensions.")
#     end

#     if dim > 1
#         Vb_copy = copy(val_fnct.Vb)

#         # permute dimensions
#         permute_vec = copy(val_fnct.state_space_length)
#         permute_vec[1] = dim
#         permute_vec[dim] = 1
#         Vb_copy = permutedims(Vb_copy, permute_vec)
#     else
#         Vb_copy = val_fnct.Vb
#     end
#     cds = size(Vb_copy, 1) # cds for current dimension size
#     for i = 1:prod(grid.state_space_length[1:end .!= dim])
#         Vb_copy[1 + (i-1)*cds] = bndy[i]
#     end
#     if dim > 1
#         val_fnct.Vb = permutedims(Vb_copy, permute_vec)
#     end
# end
