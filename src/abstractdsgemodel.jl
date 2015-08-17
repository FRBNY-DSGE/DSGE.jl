abstract AbstractDSGEModel

Base.getindex(m::AbstractDSGEModel, i::Integer) =
	(j = i - length(m.parameters)) > 0 ? m.steady_state[j] : m.parameters[i]
Base.getindex(m::AbstractDSGEModel, k::Symbol) =
	(j = (i = m.keys[k]) - length(m.parameters)) > 0 ? m.steady_state[j] : m.parameters[i]

Base.setindex!(m::AbstractDSGEModel, value, i::Integer) =
	(j = i - length(m.parameters)) > 0 ? setindex!(m.steady_state, value, j) : setindex!(m.parameters, value, i)
Base.setindex!(m::AbstractDSGEModel, value, k::Symbol) =
	(j = (i = m.keys[k]) - length(m.parameters)) > 0 ? (@inbounds m.steady_state[j] = value) : (@inbounds m.parameters[i] = value)

function Base.show{T<:AbstractDSGEModel}(io::IO, m::T)
    @printf io "%s" "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "%s %s\n" "model: " T
    @printf io "%s             %i\n" "no. states:" states(m)
    @printf io "%s %i\n" "no. anticipated shocks:" anticipated_shocks(m)
    @printf io "%s   %i\n" "no. anticipated lags:" anticipated_lags(m)
    @printf io "%s\n %s\n" "description:" description(m)
end

# Number of anticipated policy shocks
anticipated_shocks(m::AbstractDSGEModel) = m.anticipated_shocks

# Padding for nant
anticipated_shocks_padding(m::AbstractDSGEModel) = m.anticipated_shocks_padding

# Number of periods back we should start incorporating zero bound expectations
# ZLB expectations should begin in 2008 Q4
anticipated_lags(m::AbstractDSGEModel) = m.anticipated_lags

# TODO: This should be set when the data are read in
# Number of presample periods
presample_periods(m::AbstractDSGEModel) = m.presample_periods

# Number of a few things that are useful apparently
states(m::AbstractDSGEModel)      = length(m.endogenous_states)
augmented_states(m::AbstractDSGEModel)  = states(m) + length(m.endogenous_states_postgensys)
exogenous_shocks(m::AbstractDSGEModel)   = length(m.exogenous_shocks)
expected_shocks(m::AbstractDSGEModel)   = length(m.expected_shocks)
equilibrium_conditions(m::AbstractDSGEModel)     = length(m.equilibrium_conditions)
observables(m::AbstractDSGEModel) = length(m.observables)
