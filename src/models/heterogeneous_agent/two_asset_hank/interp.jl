"""
```
function interp(x::Vector{T}, x_knot::Vector{T}) where T<:AbstractFloat
```
Compute linear interpolation over 1D grid. It extrapolates outside of
knot points. (SeHyoun Ahn, March 2017)

# Arguments:
- `x_fine`: x grid to interpolate to
- `x_knot`: x knot points to interpolate from

# Output:
`V`: matrix giving interpolation values

# Example:
```
x_fine = linspace(-1,1,300)'
x_knot = linspace(0,2,20)'
V = interp(x_fine, x_knot)
plot(V * (x_knot.^5 + x_knot * 0.5))
```
"""
function interp(x::Vector{T}, x_knot::Vector{T}) where {T<:AbstractFloat}
    loc = sum(broadcast(-, x, x_knot') .>= 0, dims=2)
    loc = max.(min.(loc, length(x_knot) - 1), 1)

    t = (x - x_knot[loc]) ./ (x_knot[loc+1] - x_knot[loc])
    ind_x = 1:length(x)

    i = vec(repeat(ind_x,2,1))
    j = [loc; loc + 1]
    v = [(1 - t); t]

    return SparseArrays.sparse(i, j, v, length(x), length(x_knot))
end

"""
```
function interp(x_fine::Vector{T}, y_fine::Vector{T},
                x_knot::Vector{T}, y_knot::Vector{T}) where {T<:Float64}
```
Compute linear interpolation over 2D grid. It extrapolates outside of
knot points. (SeHyoun Ahn, March 2017)

# Arguments
 - `x_fine::Vector{T}`: x grid to interpolate to
 - `y_fine::Vector{T}`: y grid to interpolate to
 - `x_knot::Vector{T}`: x knot points to interpolate from
 - `y_knot::Vector{T}`: y knot points to interpolate from

# Output
 - `V`: sparse matrix giving interpolation values

# Example
```
x_fine = linspace(0,2,80)'
y_fine = linspace(-1,1,100)'
x_knot = linspace(0,1,10)'
y_knot = linspace(0,1,10)'

V = interp(x_fine, y_fine, x_knot, y_knot)

z = broadcast(+, x_knot.^3 + exp(-x_knot), (y_knot' - 0.5).^2)
surf(x_fine,y_fine,reshape(V * vec(z), 80, 100)')
```
"""
function interp(x_fine::Vector{T}, y_fine::Vector{T},
                x_knot::Vector{T}, y_knot::Vector{T}) where {T<:Float64}

    n_x      = length(x_knot)
    ind_fine = 1:length(x_fine) * length(y_fine)

    loc_x = sum(broadcast(-, x_fine, x_knot') .>= 0, dims=2)
    loc_x = vec(max.(min.(loc_x, n_x - 1), 1))

    loc_y = sum(broadcast(-, y_fine, y_knot') .>= 0, dims=2)
    loc_y = max.(min.(loc_y, length(y_knot) - 1), 1)'

    t_x = (x_fine .- x_knot[loc_x]) ./ (x_knot[loc_x .+ 1] .- x_knot[loc_x])
    t_y = ((y_fine .- vec(y_knot[loc_y])) ./
           (vec(y_knot[loc_y .+ 1]) .- vec(y_knot[loc_y])))'

    j_SW = broadcast(+, loc_x,      n_x.*(loc_y .- 1))
    j_SE = broadcast(+, loc_x .+ 1, n_x.*(loc_y .- 1))
    j_NW = broadcast(+, loc_x,      n_x.*(loc_y))
    j_NE = broadcast(+, loc_x .+ 1, n_x.*(loc_y))
    v_SW = broadcast(*, (1 .- t_x), (1 .- t_y))
    v_SE = broadcast(*, (t_x),      (1 .- t_y))
    v_NW = broadcast(*, (1 .- t_x), (t_y))
    v_NE = broadcast(*, (t_x),      (t_y))

    i = vec(repeat(ind_fine, 4, 1))
    j = [vec(j_SW); vec(j_SE); vec(j_NW); vec(j_NE)]
    v = [vec(v_SW); vec(v_SE); vec(v_NW); vec(v_NE)]

    return SparseArrays.sparse(i, j, v, length(x_fine) * length(y_fine), n_x * length(y_knot))
end
