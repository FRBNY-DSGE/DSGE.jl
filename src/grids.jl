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
    Grid(reshape(grid, size(grid)[1]), weights, scale)
end

####################
# Quadrature rules
####################
# Enforce that quadrature rules with additional keyword arguments get mapped to
# A pre-populated version of that rule that only accepts 3 arguments:
# Lower bound, upper bound, and number of points
function uniform_quadrature(lower_bound::T, upper_bound::T, n_points::Int;
                            scale::T = 1) where {T<:Real}
    grid = collect(range(lower_bound, stop = upper_bound, length = n_points))
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
    return grid', dropdims(weights', dims=2)
end

function curtis_clenshaw_quadrature(kind::Int = 2)
    return (lb, ub, n) -> curtis_clenshaw_quadrature(lb, ub, n, kind = kind)
end

# Tauchen86 developed a method for "finding a discrete-valued Markov chain
# whose sample paths approximate well those of a vector autoregression"
# Hence the method below is the Tauchen86 method applied to a general AR(1).
function tauchen86(μ::AbstractFloat,ρ::AbstractFloat,σ::AbstractFloat,n::Int64,λ::AbstractFloat)
    #output is xgrid, xprob
    # x_t+1 = μ + ρ x_t + σ e_{t+1}, e_t+1 ∼ N(0,1)
    xhi = μ/(1-ρ) + λ*sqrt(σ^2/(1-ρ)^2)
    xlo = μ/(1-ρ) - λ*sqrt(σ^2/(1-ρ)^2)
    xgrid = zeros(n);
    xscale =(xhi-xlo)/(n-1)

    for i=1:n
        xgrid[i] = xlo + xscale*(i-1)
    end
    m=zeros(n-1);
    for i=1:n-1
        m[i] = (xgrid[i]+xgrid[i+1])/2
    end
    xprob = zeros(n,n) # xprob[i,j] = Pr(x_t+1=xgrid[j]|x_t=xgrid[i])
    for i=1:n # this is the state today
        normie = Normal(μ+ρ*xgrid[i],σ)
        normpdf(x) = pdf.(normie,x)
        normcdf(x) = cdf.(normie,x)
        for j=2:n-1
            xprob[i,j] = normcdf(m[j]) - normcdf(m[j-1])
        end
        xprob[i,1] = normcdf(m[1])
        xprob[i,n] = 1 - normcdf(m[n-1])
    end
    xprob = xprob./sum(xprob, dims = 2) # make sure the rows sum to 1
    return ( xgrid,xprob, xscale )
end

# However, Tauchen86 can also be applied to "degenerate" AR(1)s where ρ = 0.
# i.e. i.i.d noise with mean μ, std. σ, but the way the returned
# Markov transition matrix, xprob, is used in the code, assumes it to be
# 1-dimensional, since the second-dimension of the Markov transition matrix is trivial
# (since the process is i.i.d).
function tauchen86(μ::AbstractFloat, σ::AbstractFloat, n::Int64, λ::AbstractFloat)
    xgrid, xprob, xscale = tauchen86(μ, 0., σ, n, λ)
    return xgrid, vec(xprob[:, 1]), xscale
end

####################
# Grid-based utils
####################
function get_grid(m::AbstractDSGEModel, grid_name::Symbol)
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
