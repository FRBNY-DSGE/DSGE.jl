# Adds the unprimed model states to the dictionary of ranges to properly
# fill the Jacobian matrix
# Note: This is not augmentation in the usual sense, which incorporates lags of
# model state variables after the transition equation has been solved for
# but rather, augmenting the model states with lags prior to solution since the
# solution method requires lags
function augment_model_states(endo::AbstractDict{Symbol, UnitRange}, n_model_states::Int64)
    endo_aug = deepcopy(endo)
    for (state::Symbol, inds::UnitRange) in endo
        unprimed_state = unprime(state)
        unprimed_inds::UnitRange  = reindex_unprimed_model_states(inds, n_model_states)
        endo_aug[unprimed_state] = unprimed_inds
    end

    # Ensure all ranges are consecutive
    endo_ranges = endo_aug.vals
  #=  for i in 1:length(endo_ranges)
        if i == 1
            @assert endo_ranges[1].start == 1
        else
            @assert endo_ranges[i].start == endo_ranges[i-1].stop + 1
        end
    end =#

    return endo_aug
end

# For the RANK equivalent of the HANK models
function augment_model_states(endo::AbstractDict{Symbol, Int64}, n_model_states::Int64)
    endo_aug = deepcopy(endo)
    for (state::Symbol, ind::Int) in endo
        unprimed_state           = unprime(state)
        endo_aug[unprimed_state] = ind + n_model_states
    end

    return endo_aug
end

# The reason for this function is that the canonical form for the Klein solution method is
# as follows:
# E_t A([x_{t+1}, y_{t+1}]) = B[x_t, y_t]
# Where we only need to track [x_{t+1}, y_{t+1}] as states when
# we transform the system from canonical form to state space form
# Hence because we only keep track of the t+1 indexed variables, which we denote
# with a ′, we need a way of calculating the indices in the Jacobian corresponding
# to the t indexed variables.
function reindex_unprimed_model_states(inds::UnitRange, n_model_states::Int64)
    return inds .+ n_model_states
end

# Returns a model state variable symbol without the prime
function unprime(state::Symbol)
    return Symbol(replace(string(state), "′" => ""))
end

######################################################################################
# For normalizing the states
function normalize_model_state_indices!(m::AbstractModel)
    normalized_model_states = m.normalized_model_states

    endo                     = m.endogenous_states
    model_state_keys         = collect(keys(endo))
    normalized_model_states = intersect(model_state_keys, normalized_model_states) #findall((in)(normalized_model_states), model_state_keys)

    state_inds                  = get_setting(m, :state_indices)
    jump_inds                   = get_setting(m, :jump_indices)
    state_normalization_factor  = get_setting(m, :backward_looking_states_normalization_factor)
    jump_normalization_factor   = get_setting(m, :jumps_normalization_factor)

    m.endogenous_states = normalize(m, endo, model_state_keys, normalized_model_states,
                                    state_inds, jump_inds,
                                    state_normalization_factor, jump_normalization_factor)
end

# Shift a UnitRange type down by an increment
# If this UnitRange is the first_range in the group being shifted, then do not
# subtract the increment from the start
# e.g.
# If you have 1:80 as your first range, then shift(1:80, -1; first_range = true)
# should return 1:79
# However, if you have 81:160 as range, that is not your first range, then the
# function should return 80:159
function shift(inds::UnitRange, increment::Int64; first_range::Bool = false)
    if first_range
        return UnitRange(inds.start, inds.stop + increment)
    else
        return UnitRange(inds.start + increment, inds.stop + increment)
    end
end

# normalize the distributional states located in indices specified by normalized_state_inds
function normalize(m::AbstractModel,
                   endo::AbstractDict{Symbol, UnitRange},
                   model_state_keys::Vector{Symbol},
                   normalized_model_states::Vector{Symbol},
                   state_inds::AbstractArray{Int64},
                   jump_inds::AbstractArray{Int64},
                   state_normalization_factor::Int64,
                   jump_normalization_factor::Int64)
    n_model_state_vars = length(endo)

    # For each model state that needs to be normalized...
    for state in normalized_model_states
        normalization_factor =
        if state in get_setting(m, :states)
            state_normalization_factor
        elseif state in get_setting(m, :jumps)
            jump_normalization_factor
        end
        # Subtract 1 from both the beginning and end of the UnitRange...
        for state in model_state_keys
            inds = endo[state]
            # except for first UnitRange in a given normalization
            first_range = (state==get_setting(m, :states)[1])
            endo[state] = shift(inds, -normalization_factor, first_range = first_range)
        end
        #=for j in i:n_model_state_vars
            inds = endo[model_state_keys[j]]
            # Except for the first UnitRange in a given normalization
            first_range = (j == i)
            endo[model_state_keys[j]] = shift(inds, -normalization_factor, first_range = first_range)
        end=#
    end

    # Ensure all ranges are consecutive
    endo_ranges = endo.vals
   #= for i in 1:n_model_state_vars
        if i == 1
            @assert endo_ranges[1].start == 1
        else
            @assert endo_ranges[i].start == endo_ranges[i-1].stop + 1
        end
    end =#

    return endo
end


"""
```
stack_indices(key_dict::OrderedDict, keys::Vector{Symbol})
```
Stacks the indices of a OrderedDict{Symbol, AbstractRange} into a Vector of Ints defined by the Ranges in the dictionary.
"""
function stack_indices(key_dict::AbstractDict, keys::Vector{Symbol})
    indices = Vector{Int64}(undef, 0)
    for i in getindex.(Ref(key_dict), keys)
        indices = [indices; i]
    end
    return indices
end
