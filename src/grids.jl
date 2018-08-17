###############################
# Type Definition/Constructors
###############################
type Grid
    points::Vector{Float64}
    weights::Vector{Float64}
    scale::Float64
    Grid(points, weights, scale) = sum(weights) ≈ scale ? new(points, weights, scale) : error("scaled weights do not sum up properly")
end

# Constructor utilizing a custom weight calculation function
function Grid{T<:Real}(quadrature::Function,
                       lower_bound::T, upper_bound::T,
                       n_points::Int; scale::Float64 = 1.)
    grid, weights = quadrature(lower_bound, upper_bound, n_points)
    Grid(grid, weights, scale)
end

####################
# Quadrature rules
####################
# Enforce that quadrature rules with additional keyword arguments get mapped to
# A pre-populated version of that rule that only accepts 3 arguments:
# Lower bound, upper bound, and number of points
function uniform_quadrature{T<:Real}(lower_bound::T, upper_bound::T, n_points::Int;
                                     scale::T = 1)
    grid = collect(linspace(lower_bound, upper_bound, n_points))
    weights = fill(scale/n_points, n_points)
    return grid, weights
end

function uniform_quadrature{T<:Real}(scale::T = 1)
    return (lb, ub, n) -> uniform_quadrature(lb, ub, n, scale = scale)
end

# NOTE: Should modify chebpts function to return cleaner types...
# Kind indicates the "kind" of the chebyshev polynomial grid
# The terminology is generally, "Chebyshev polynomials of the Nth kind"
function curtis_clenshaw_quadrature{T<:Real}(lower_bound::T, upper_bound::T,
                                             n_points::Int; kind::Int = 2)
    grid, weights = chebpts(n_points, lower_bound, upper_bound, kind)
    return grid', squeeze(weights', 2)
end

function curtis_clenshaw_quadrature(kind::Int = 2)
    return (lb, ub, n) -> curtis_clenshaw_quadrature(lb, ub, n, kind = kind)
end

####################
# Grid-based utils
####################
function get_grid(m::AbstractModel, grid_name::Symbol)
    return m.grids[grid_name]
end

function quadrature_sum{T<:Real}(x::Vector{T}, grid::Grid)
    return sum(grid.weights .* x .* grid.points)
end

function quadrature_sum{T<:Real}(grid::Grid, x::Vector{T})
    return sum(x, grid)
end

# # Defining arithmetic on/standard function evaluation of grids
# for op in (:(Base.:+),
           # :(Base.:-),
           # :(Base.:*),
           # :(Base.:/))
           # # Exponentiation also omitted because of broadcasting issues

    # @eval ($op)(g::Grid, x::Integer)           = ($op)(g.points, x)
    # @eval ($op)(g::Grid, x::Number)            = ($op)(g.points, x)
    # if op == :(Base.:/)
        # op = :(Base.:./) # Will not work because the syntax for broadcasting has been
                         # # updated to broadcast(/, ...)?
        # @eval ($op)(x::Integer, g::Grid)       = ($op)(x, g.points)
        # @eval ($op)(x::Number, g::Grid)        = ($op)(x, g.points)
    # end
# end

# for f in (:(Base.:-),
          # :(Base.:<),
          # :(Base.:>),
          # :(Base.:<=),
          # :(Base.:>=))

    # @eval ($f)(g::Grid) = ($f)(g.points)

    # # Similar issues with broadcasting...
    # if f != :(Base.:-)
        # @eval ($f)(g::Grid, x::Number) = ($f)(g.points, x)
    # end
# end
