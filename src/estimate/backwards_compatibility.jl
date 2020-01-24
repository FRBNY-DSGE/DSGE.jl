"""
```
function old_to_new_cloud(cloud::DSGE.Cloud)::SMC.Cloud
```
Converter between types.
"""
function old_to_new_cloud(cloud::DSGE.Cloud)::SMC.Cloud
    return SMC.Cloud(cloud.particles, cloud.tempering_schedule,
                     cloud.ESS, cloud.stage_index, cloud.n_Φ, cloud.resamples,
                     cloud.c, cloud.accept, cloud.total_sampling_time)
end

"""
```
function Cloud(c::ParticleCloud)
function Cloud(c::Cloud)
```
Returns a Cloud type.
"""
function SMC.Cloud(cloud::ParticleCloud)
    return SMC.Cloud(hcat(Matrix{Float64}(get_vals(cloud)'), get_loglh(cloud), get_logprior(cloud),
                      get_old_loglh(cloud), get_accept(cloud), get_weights(cloud)),
                 cloud.tempering_schedule, cloud.ESS, cloud.stage_index, cloud.n_Φ,
                 cloud.resamples, cloud.c, cloud.accept, cloud.total_sampling_time)
end
function SMC.Cloud(cloud::SMC.Cloud)
    return cloud
end
