abstract AbstractDSGEModel{T<:FloatingPoint}

function Base.show{T<:AbstractDSGEModel}(io::IO, m::T)
    @printf io "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "%s\n" T
    @printf io "no. states:             %i\n" num_states(m)
    @printf io "no. anticipated shocks: %i\n" num_anticipated_shocks(m)
    @printf io "no. anticipated lags:   %i\n" num_anticipated_lags(m)
    @printf io "description:\n %s\n"          description(m)
end

# TODO consider stacking all parameters in a single vector. Alternately, all fixed
# parameters can be added to the normal parameters vector at a (potentially negligible)
# performance hit.
@inline function Base.getindex(m::AbstractDSGEModel, i::Integer)
    if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]q
    end
end

# need to define like this so we can disable bounds checking
@inline function Base.getindex(m::AbstractDSGEModel, k::Symbol)
    i = m.keys[k]
    @inbounds if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]
    end
end

@inline function Base.setindex!(m::AbstractDSGEModel, value, i::Integer)
    if i <= (j = length(m.parameters))
        return setindex!(m.parameters, value, i)
    else
        return setindex!(m.steady_state, value, i-j)
    end
end
Base.setindex!(m::AbstractDSGEModel, value, k::Symbol) = Base.setindex!(m, value, m.keys[k])

# syntax for adding a prameter to a model m <= parameter
function (<=){T}(m::AbstractDSGEModel{T}, p::AbstractParameter{T})
    @assert !in(p.key, m.keys) "Key $(p.key) is already present in DSGE model"

    new_param_index = length(m.keys) + 1

    # grow parameters and add the parameter
    push!(m.parameters, p)

    # add parameter location to dict
    setindex!(m.keys, new_param_index, p.key)
end

function (<=){T}(m::AbstractDSGEModel{T}, k::Symbol)
    @assert !in(k, m.keys) "Key $(k) is already present in DSGE model"

    new_param_index = length(m.keys) + 1

    # grow steady_state values with a zero
    push!(m.steady_state, zero(T))

    # add parameter location to dict
    setindex!(m.keys, new_param_index, k)
end

function (<=)(m::AbstractDSGEModel, vec::Vector{Symbol})
    for k in vec
        m <= k
    end
end

Distributions.logpdf(m::AbstractDSGEModel) = logpdf(m.parameters)
Distributions.pdf(m::AbstractDSGEModel) = exp(logpdf(m))

# Number of anticipated policy shocks
num_anticipated_shocks(m::AbstractDSGEModel) = m.num_anticipated_shocks

# Padding for nant
num_anticipated_shocks_padding(m::AbstractDSGEModel) = m.num_anticipated_shocks_padding

# Number of periods back we should start incorporating zero bound expectations
# ZLB expectations should begin in 2008 Q4
num_anticipated_lags(m::AbstractDSGEModel) = m.num_anticipated_lags

# TODO: This should be set when the data are read in
# Number of presample periods
num_presample_periods(m::AbstractDSGEModel) = m.num_presample_periods

# Number of a few things that are useful apparently
num_states(m::AbstractDSGEModel)                 = length(m.endogenous_states)
num_states_augmented(m::AbstractDSGEModel)       = num_states(m) + length(m.endogenous_states_postgensys)
num_shocks_exogenous(m::AbstractDSGEModel)       = length(m.exogenous_shocks)
num_shocks_expectational(m::AbstractDSGEModel)   = length(m.expected_shocks)
num_equilibrium_conditions(m::AbstractDSGEModel) = length(m.equilibrium_conditions)
num_observables(m::AbstractDSGEModel)            = length(m.observables)
num_parameters(m::AbstractDSGEModel)             = length(m.parameters)
num_parameters_fixed(m::AbstractDSGEModel)       = length(m.parameters_fixed)
num_parameters_steady_state(m::AbstractDSGEModel)= length(m.steady_state)
num_parameters_free(m::AbstractDSGEModel)        = sum([!α.fixed for α in m.parameters])

# Paths to where input/output/results data are stored
savepath(m::AbstractDSGEModel)  = normpath(m.savepath)
inpath(m::AbstractDSGEModel)    = normpath(joinpath(m.savepath, "input_data/"))
outpath(m::AbstractDSGEModel)   = normpath(joinpath(m.savepath, "output_data/"))
tablepath(m::AbstractDSGEModel) = normpath(joinpath(m.savepath, "results/tables/"))
plotpath(m::AbstractDSGEModel)  = normpath(joinpath(m.savepath, "results/plots/"))
logpath(m::AbstractDSGEModel)   = normpath(joinpath(m.savepath, "logs/"))

# TODO is there a better place for these? They do depend on AbstractDSGEModel type.
function tomodel!{T<:FloatingPoint}(m::AbstractDSGEModel, values::Vector{T})
    tomodel!(values, m.parameters)
    return steadystate!(m)
end
function update!{T<:FloatingPoint}(m::AbstractDSGEModel, values::Vector{T})
    return update!(m.parameters, values)
end

