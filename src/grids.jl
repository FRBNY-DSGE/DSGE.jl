import Distributions: Normal, pdf, cdf

###############################
# Type Definition/Constructors
###############################
mutable struct Grid
    points::Vector{Float64}
    weights::Vector{Float64}
    scale::Float64
    Grid(points, weights, scale) = sum(weights) ≈ scale ? new(points, weights, scale) : error("scaled weights do not sum up properly")
end

# Constructor utilizing a custom weight calculation function
function Grid(quadrature::Function,
              lower_bound::T, upper_bound::T,
              n_points::Int; scale::Float64 = 1.) where {T<:Real}
    grid, weights = quadrature(lower_bound, upper_bound, n_points)
    Grid(grid, weights, scale)
end

####################
# Quadrature rules
####################
# Enforce that quadrature rules with additional keyword arguments get mapped to
# A pre-populated version of that rule that only accepts 3 arguments:
# Lower bound, upper bound, and number of points
function uniform_quadrature(lower_bound::T, upper_bound::T, n_points::Int;
                            scale::T = 1) where {T<:Real}
    grid = collect(range(lower_bound, upper_bound, n_points))
    weights = fill(scale/n_points, n_points)
    return grid, weights
end

function uniform_quadrature(scale::T = 1) where {T<:Real}
    return (lb, ub, n) -> uniform_quadrature(lb, ub, n, scale = scale)
end

# NOTE: Should modify chebpts function to return cleaner types...
# Kind indicates the "kind" of the chebyshev polynomial grid
# The terminology is generally, "Chebyshev polynomials of the Nth kind"
function curtis_clenshaw_quadrature(lower_bound::T, upper_bound::T,
                                    n_points::Int; kind::Int = 2) where {T<:Real}
    grid, weights = chebpts(n_points, lower_bound, upper_bound, kind)
    return grid', squeeze(weights', 2)
end

function curtis_clenshaw_quadrature(kind::Int = 2)
    return (lb, ub, n) -> curtis_clenshaw_quadrature(lb, ub, n, kind = kind)
end

function tauchen86(μ::AbstractFloat,σ::AbstractFloat,n::Int64,λ::AbstractFloat)
    # output is xgrid, xprob
    xhi = μ + λ*σ;
    xlo = μ - λ*σ;
    xgrid = zeros(n);
    xscale =(xhi-xlo)/(n-1)

    for i=1:n
        xgrid[i] = xlo + xscale*(i-1);
    end
    m = zeros(n-1);
    for i=1:n-1
        m[i] = (xgrid[i]+xgrid[i+1])/2;
    end
    xprob = zeros(n);
    normie = Normal(μ,σ)
    normpdf(x) = pdf.(normie,x)
    normcdf(x) = cdf.(normie,x)
    for j=2:n-1
        xprob[j] = normcdf(m[j]) - normcdf(m[j-1]);
    end
    xprob[1] = normcdf(m[1]);
    xprob[n] = 1 - normcdf(m[n-1]);
    return ( xgrid,xprob, xscale )
end

####################
# Grid-based utils
####################
function get_grid(m::AbstractModel, grid_name::Symbol)
    return m.grids[grid_name]
end

function quadrature_sum(x::Vector{T}, grid::Grid) where {T<:Real}
    return sum(grid.weights .* x .* grid.points)
end

function quadrature_sum(grid::Grid, x::Vector{T}) where {T<:Real}
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
