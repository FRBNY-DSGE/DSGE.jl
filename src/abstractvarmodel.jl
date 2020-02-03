"""
```
AbstractVARModel{T} <: AbstractModel{T}
```

The AbstractVARModel is defined as a subtype of AbstractModel to implement
vector autoregression methods.
"""
abstract type AbstractVARModel{T} <: ModelConstructors.AbstractModel{T} end

function Base.show(io::IO, m::AbstractVARModel)
    @printf io "VAR Model\n"
    @printf io "observables:   %i\n" observables(m)
    @printf io "data_vintage:  %i\n" data_vintage(m)
    @printf io "description:\n %s\n" description(m)
end

# """
# ```
# AbstractDSGEVARModel{T} <: AbstractVARModel{T}
# ```

# The AbstractDSGEVARModel is defined as a subtype of AbstractVARModel to implement
# DSGE-VAR methods.
# """
# abstract type AbstractDSGEVARModel{T} <: AbstractVARModel{T} end

# function Base.show(io::IO, m::AbstractDSGEVAR)
#     @printf io "DSGE-VAR Model\n"
#     @printf io "observables:      %i\n" observables(m)
#     @printf io "data_vintage:     %i\n" data_vintage(m)
#     @printf io "DSGE model:       %s\n" spec(get_dsge(m))
#     @printf io "DSGE description: %s\n" description(get_dsge(m))
# end

# get_dsge(m::AbstractDSGEVAR) = m.dsge
