"""
```
initial_draw(m::AbstractModel, data::Array{Float64}, c::ParticleCloud)
```

Draw from a general starting distribution (set by default to be from the prior) to initialize the SMC algorithm.
Returns a tuple (logpost, loglh) and modifies the particle objects in the particle cloud in place.

"""
function initial_draw(m::AbstractModel, data::Array{Float64}, c::ParticleCloud)
    dist_type = get_setting(m, :initial_draw_source)
    if dist_type == :normal
        params = zeros(n_parameters(m))
        hessian = zeros(n_parameters(m),n_parameters(m))
        try
            file = h5open(rawpath(m, "estimate", "paramsmode.h5"), "r")
            params = read(file,"params")
            close(file)
        catch
            throw("There does not exist a valid mode file at "*rawpath(m,"estimate","paramsmode.h5"))
        end

        try
            file = h5open(rawpath(m, "estimate", "hessian.h5"), "r")
            hessian = read(file, "hessian")
            close(file)
        catch
            throw("There does not exist a valid hessian file at "*rawpath(m,"estimate","hessian.h5"))
        end

        S_diag, U = eig(hessian)
        big_eig_vals = find(x -> x>1e-6, S_diag)
        rank = length(big_eig_vals)
        n = length(params)
        S_inv = zeros(n, n)
        for i = (n-rank+1):n
            S_inv[i, i] = 1/S_diag[i]
        end
        hessian_inv = U*sqrt(S_inv)
        dist = DSGE.DegenerateMvNormal(params, hessian_inv)
    end

    n_part = length(c)
    draws =
    dist_type == :prior ? rand(m.parameters, n_part) : rand(dist, n_part)

    loglh = zeros(n_part)
    logpost = zeros(n_part)
    for i in 1:size(draws)[2]
        success = false
        while !success
            try
                update!(m, draws[:, i])
                loglh[i] = likelihood(m, data)
                logpost[i] = prior(m)
            catch
                draws[:, i] =
                dist_type == :prior ? rand(m.parameters, 1) : rand(dist, 1)
                continue
            end
            success = true
        end
    end
    update_draws!(c, draws)
    update_loglh!(c, loglh)
    update_logpost!(c, logpost)
end

"""
```
mvnormal_mixture_draw{T<:AbstractFloat}(θ_old, σ; cc, α, θ_prop)
```

Create a `DegenerateMvNormal` distribution object, `d`, from a parameter vector, `p`, and a
standard deviation matrix (obtained from SVD), `σ`.

Generate a draw from the mixture distribution of:
1. A `DegenerateMvNormal` centered at θ_old with the standard deviation matrix `σ`, scaled by `cc^2` and with mixture proportion `α`.
2. A `DegenerateMvNormal` centered at the same mean, but with a standard deviation matrix of the diagonal entries of `σ` scaled by `cc^2` with mixture proportion `(1 - α)/2`.
3. A `DegenerateMvNormal`  with the same standard deviation matrix `σ` but centered at the new proposed mean, `θ_prop`, scaled by `cc^2`, and with mixture proportion `(1 - α)/2`.

If no `θ_prop` is given, but an `α` is specified, then the mixture will consist of `α` of
the standard distribution and `(1 - α)` of the diagonalized distribution.

### Arguments
- `θ_old::Vector{T}`: The mean of the desired distribution
- `σ::Matrix{T}`: The standard deviation matrix of the desired distribution

### Keyword Arguments
- `cc::T`: The standard deviation matrix scaling factor
- `α::T`: The mixing proportion
- `θ_prop::Vector{T}`: The proposed parameter vector to be used as part of the mixture distribution, set by default to be the weighted mean of the particles, prior to mutation.

### Outputs
- `θ_new::Vector{T}`: The draw from the mixture distribution to be used as the MH proposed step
- `new mixture_density::T`: The mixture density conditional on θ_old evaluated at `θ_new` to be used in calculating the MH move probability
- `old mixture_density::T`: The mixture density conditional on θ_new evaluated at `θ_old` to be used in calculating the MH move probability

"""
function mvnormal_mixture_draw{T<:AbstractFloat}(θ_old::Vector{T}, σ::Matrix{T};
                                                 cc::T = 1.0, α::T = 1.,
                                                 θ_prop::Vector{T} = θ_old)
    @assert 0 <= α <= 1

    # Create mixture distribution conditional on the previous parameter value, θ_old
    d_old = DegenerateMvNormal(θ_old, cc*σ)
    d_diag_old = DegenerateMvNormal(θ_old, diagm(diag(cc*σ)))
    d_prop = θ_prop == θ_old ? d_diag : DegenerateMvNormal(θ_prop, cc*σ)
    d_mix_old = MixtureModel(DegenerateMvNormal[d_old, d_diag_old, d_prop], [α, (1 - α)/2, (1 - α)/2])

    θ_new = rand(d_mix_old)

    # Create mixture distribution conditional on the new parameter value, θ_new
    d_new = DegenerateMvNormal(θ_new, cc*σ)
    d_diag_new = DegenerateMvNormal(θ_new, diagm(diag(cc*σ)))
    d_mix_new = MixtureModel(DegenerateMvNormal[d_new, d_diag_new, d_prop], [α, (1 - α)/2, (1 - α)/2])

    # To clarify, this is not just the density of θ_new/θ_old using a given mixture
    # density function, but rather, the density of θ_new | θ_old and the density of
    # θ_old | θ_new taken with respect to their respective mixture density functions
    new_mixture_density = logpdf(d_mix_old, θ_new)
    old_mixture_density = logpdf(d_mix_new, θ_old)

    return θ_new, new_mixture_density, old_mixture_density
end

function compute_ESS{T<:AbstractFloat}(loglh::Vector{T}, current_weights::Vector{T},
                                       ϕ_n::T, ϕ_n1::T; use_CESS::Bool = false)
    incremental_weights = exp.((ϕ_n - ϕ_n1)*loglh)
    new_weights = current_weights.*incremental_weights
    if use_CESS
        N   = length(incremental_weights)
        ESS = N*sum(current_weights .* incremental_weights)^2/sum(current_weights .*
                                                                  incremental_weights.^2)
    else
        normalized_weights = new_weights/sum(new_weights)
        ESS = 1/sum(normalized_weights.^2)
    end
    return ESS
end

function init_stage_print(cloud::ParticleCloud;
                          verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    if use_fixed_schedule
        println("--------------------------")
            println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
    else
        println("--------------------------")
            println("Iteration = $(cloud.stage_index)")
    end
	println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
	println("--------------------------")
        println("c = $(cloud.c)")
	println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
	    end
    end
end

function end_stage_print(cloud::ParticleCloud;
                         verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    total_sampling_time_minutes = cloud.total_sampling_time/60
    if use_fixed_schedule
        expected_time_remaining_sec = (cloud.total_sampling_time/cloud.stage_index)*(cloud.n_Φ - cloud.stage_index)
        expected_time_remaining_minutes = expected_time_remaining_sec/60
    end

    println("--------------------------")
    if use_fixed_schedule
        println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
        println("time elapsed: $(round(total_sampling_time_minutes, 4)) minutes")
        println("estimated time remaining: $(round(expected_time_remaining_minutes, 4)) minutes")
    else
        println("Iteration = $(cloud.stage_index)")
        println("time elapsed: $(round(total_sampling_time_minutes, 4)) minutes")
    end
    println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
    println("--------------------------")
        println("c = $(cloud.c)")
        println("accept = $(cloud.accept)")
        println("ESS = $(cloud.ESS[cloud.stage_index])   ($(cloud.resamples) total resamples.)")
    println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
        end
    end
end
