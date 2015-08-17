abstract AbstractDSGEModel

Base.getindex(m::AbstractDSGEModel, i::Integer) = m.par[i]
Base.getindex(m::AbstractDSGEModel, keyword::Symbol) = m.par[m.parkeys[keyword]]

Base.setindex!(m::AbstractDSGEModel, value, i::Integer) = setindex!(m.par, value, i)
Base.setindex!(m::AbstractDSGEModel, value, keyword::Symbol) = setindex!(m.par, value, m.parkeys[keyword])

function Base.show{T<:AbstractDSGEModel}(io::IO, m::T)
    @printf io "%s" "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "%s %s\n" "model: " T
    @printf io "%s             %i\n" "# states:" states(m)
    @printf io "%s %i\n" "# anticipated shocks:" anticipated_shocks(m)
    @printf io "%s   %i\n" "# anticipated lags:" anticipated_lags(m)
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
