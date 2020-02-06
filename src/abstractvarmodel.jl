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

get_parameters(m::AbstractVARModel) = m.parameters

"""
```
AbstractDSGEVARModel{T} <: AbstractVARModel{T}
```
The AbstractDSGEVARModel is defined as a subtype of AbstractVARModel to implement
DSGE-VAR methods.
"""
abstract type AbstractDSGEVARModel{T} <: AbstractVARModel{T} end

function Base.show(io::IO, m::AbstractDSGEVARModel)
    @printf io "DSGE-VAR Model\n"
    @printf io "observables:      %i\n" n_observables(m)
    @printf io "data_vintage:     %s\n" data_vintage(m.dsge)
    @printf io "DSGE model:       %s\n" spec(get_dsge(m))
    @printf io "DSGE description: %s\n" description(get_dsge(m))
end

# Interface for accessing parametesr
get_parameters(m::AbstractDSGEVARModel) = m.dsge.parameters


# Forecast input files
function forecast_input_file_overrides(m::AbstractDSGEVARModel)
    return forecast_input_file_overrides(m.dsge)
end

function get_forecast_input_file(m::AbstractDSGEVARModel)
    return get_forecast_input_file(m.dsge)
end
