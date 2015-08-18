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
    @printf io "%s\n" T
    @printf io "%s             %i\n" "no. states:" num_states(m)
    @printf io "%s %i\n" "no. anticipated shocks:" num_anticipated_shocks(m)
    @printf io "%s   %i\n" "no. anticipated lags:" num_anticipated_lags(m)
    @printf io "%s\n %s\n" "description:" description(m)
end

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
num_parameters_steady_state(m::AbstractDSGEModel)= length(m.steady_state)
num_parameters_free(m::AbstractDSGEModel)        = sum([!α.fixed for α in m.parameters]) 
