abstract AbstractDSGEModel

# TODO consider stacking all parameters in a single vector. Alternately, all fixed
# parameters can be added to the normal parameters vector at a (potentially negligible)
# performance hit.
Base.getindex(m::AbstractDSGEModel, i::Integer) = begin
    if i <= num_parameters(m)
        return m.parameters[i]
    elseif i <= num_parameters(m) + num_parameters_fixed(m)
        return m.parameters_fixed[i - num_parameters(m)]
    elseif i <= num_parameters(m) + num_parameters_fixed(m) + num_parameters_steady_state(m)
        return m.steady_state[i - num_parameters(m) - num_parameters_fixed(m)]
    end
end
Base.getindex(m::AbstractDSGEModel, k::Symbol) = Base.getindex(m,m.keys[k])

Base.setindex!(m::AbstractDSGEModel, value, i::Integer) = begin 
    if i <= num_parameters(m)
        setindex!(m.parameters, value, i)
    elseif i <= num_parameters(m) + num_parameters_fixed(m)
        setindex!(m.parameters_fixed, value, i-num_parameters(m))
    elseif i <= num_parameters(m) + num_parameters_fixed(m) + num_parameters_steady_state(m)
        setindex!(m.steady_state, value, i - num_parameters(m) - num_parameters_fixed(m))
    end
end
Base.setindex!(m::AbstractDSGEModel, value, k::Symbol) = Base.setindex!(m,value,m.keys[k])

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
    update!(m.parameters, values)
    return steadystate!(m) #!!
end


## function prepareForTesting{T<:AbstractDSGEModel}(m::AbstractDSGEModel)

##     m.num_mh_simulations
##     m.num_mh_blocks
##     m.num_mh_burn
##     m.mh_thinning_step

    
## end

