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
