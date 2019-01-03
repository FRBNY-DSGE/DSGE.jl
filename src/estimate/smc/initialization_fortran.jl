"""
```
initial_draw!(m::AbstractModel, data::Matrix{Float64}, c::ParticleCloud, readin::Matrix{Float64})
```

Draw from a general starting distribution (set by default to be from the prior) to initialize the SMC algorithm.
Returns a tuple (logpost, loglh) and modifies the particle objects in the particle cloud in place.

"""
function initial_draw!(m::AbstractModel, data::Matrix{Float64}, c::ParticleCloud,
                       path::String; # RECA
                       parallel::Bool = false, use_chand_recursion::Bool = true,
                       verbose::Symbol = :low)
    # legacy - not currently necessary
    n_parts = length(c)
    loglh = zeros(n_parts)
    logpost = zeros(n_parts)

    # this is not technically the log posterior - does that matter?
    draw    = readdlm(path * "001parasim.txt")
    loglh   = readdlm(path * "001liksim.txt")
    logpost = readdlm(path * "001postsim.txt")

    update_draws!(c, draws)
    update_loglh!(c, loglh)
    update_logpost!(c, logpost)
end
