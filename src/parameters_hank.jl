"""
```
SteadyStateParameterArray{T} <: AbstractParameter{T}
```
Steady-state model parameter whose value is an Array and
 depends upon the value of other (non-steady-state)
`Parameter`s. `SteadyStateParameterArray`s must be constructed and added to an instance of a
model object `m` after all other model `Parameter`s have been defined. Once added to `m`,
`SteadyStateParameterArray`s are stored in `m.steady_state`. Their values are calculated and set
by `steadystate!(m)`, rather than being estimated directly. `SteadyStateParameterArray`s do not
require transformations from the model space to the real line or scalings for use in
equilibrium conditions.

#### Fields

- `key::Symbol`: Parameter name. Should conform to the guidelines
  established in the DSGE Style Guide.
- `value::Array{T}`: The parameter's steady-state values.
- `description::String`: Short description of the parameter's economic significance.
- `tex_label::String`: String for printing parameter name to LaTeX.
"""
type SteadyStateParameterArray{T} <: AbstractParameter{T}
    key::Symbol
    value::Array{T}
    description::String
    tex_label::String
end

"""
```
SteadyStateParameterArray{T<:Number}(key::Symbol, value::Array{T};
                                description::String = "",
                                tex_label::String = "")
```

SteadyStateParameter constructor with optional `description` and `tex_label` arguments.
"""

function SteadyStateParameterArray{T<:Number}(key::Symbol,
                                             value::Array{T};
                                             description::String = "No description available",
                                             tex_label::String = "")

    return SteadyStateParameterArray(key, value, description, tex_label)
end

# TypeError: non-boolean (BitArray{1}) used in boolean context
# gets thrown when we print the value.

function Base.show{T}(io::IO, p::SteadyStateParameterArray{T})
    @printf io "%s\n" typeof(p)
    @printf io "(:%s)\n%s\n"      p.key p.description
    @printf io "LaTeX label: %s\n"     p.tex_label
    @printf io "-----------------------------\n"
    @printf io "value:        [%+6f,...,%+6f]\n" p.value[1] p.value[end]
end

##################################
# Operator Overloading/Arithmetic
##################################

# Currently not working correctly

#=
# Defining arithmetic on/standard function evaluation of grids
for op in (:(Base.:+),
           :(Base.:-),
           :(Base.:*),
           :(Base.:/))

    @eval ($op)(g::SteadyStateParameterArray, x::Integer)        = ($op)(g.value, x)
    @eval ($op)(g::SteadyStateParameterArray, x::Number)         = ($op)(g.value, x)
    @eval ($op)(x::Integer, g::SteadyStateParameterArray)        = ($op)(x, g.value)
    @eval ($op)(x::Number, g::SteadyStateParameterArray)         = ($op)(x, g.value)
end
=#
# This should help us avoid type instability with Dual Types
#=for op in (:(Base.:+),
           :(Base.:-),
           :(Base.:*),
           :(Base.:/))

    @eval ($op)(M::Matrix{Float64}, v::Array{ForwardDiff.Dual{T,V,N},2}) where {T, V <: Float64, N} = reinterpret(ForwardDiff.Dual{T,V,N}, reshape(($op)(M, reshape(reinterpret(V,v), N+1, n))))


    @eval ($op)(M::Matrix{Float64}, v::Array{ForwardDiff.Dual{T,V,N},2}) where {T, V <: Float64, N} = reinterpret(ForwardDiff.Dual{T,V,N}, reshape(($op)(reshape(reinterpret(V,v), N+1, n), M)))
    #@eval ($op)(x::Vector, p::UnscaledOrSteadyState) = ($op)(x, p.value)
end
=#