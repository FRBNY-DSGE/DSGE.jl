# Allows for a custom thinning factor (jstep) to be specified
# If not, it pulls the jstep from the model
function thin_mh_draws(m::AbstractModel, params::Matrix{Float64}; jstep::Int64 = 1)
    jstep = jstep == 1 ? m.settings[:forecast_jstep].value : jstep
    n_total_draws, n_params = size(params)
    # Thin as usual if n_total_draws % jstep == 0
    # If it does not evenly divide, then start from the remainder+1-th index
    # and then take a thinned subset from there
    n_new_draws, offset = divrem(n_total_draws, jstep)
    params_thinned = Matrix{Float64}(n_new_draws, n_params)
    params_offset = params[offset+1:end, :]

    for (i, j) in enumerate(1:jstep:n_total_draws)
        params_thinned[i, :] = params_offset[j, :]
    end
    return params_thinned
end
