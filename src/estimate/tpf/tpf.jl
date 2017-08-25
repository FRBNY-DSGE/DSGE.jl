"""
```
tpf{S<:AbstractFloat}(m::AbstractModel, data::Array, system::System{S},
    s0::Array{S}, P0::Array{S}; verbose::Symbol=:low, include_presample::Bool=true)
```
Executes tempered particle filter.

### Inputs
- `m::AbstractModel`: model object
- `data::Array{S}`: (`n_observables` x `hist_periods`) size `Matrix{S}` of data for observables.
- `system::System{S}`: `System` object specifying state-space system matrices for model
- `s0::Array{S}`: (`n_observables` x `n_particles`) initial state vector
- `P0::Array`: (`n_observables` x `n_observavles`) initial state covariance matrix

### Keyword Arguments
- `verbose::Symbol`: indicates desired nuance of outputs. Default to `:low`.
- `include_presample::Bool`: indicates whether to include presample in periods in the returned
   outputs. Defaults to `true`.

### Outputs
- `Neff::Vector{S}`: (`hist_perdiods` x 1) vector returning inefficiency calculated per period t
- `lik::Vector{S}`: (`hist_periods` x 1) vector returning log-likelihood per period t
- `times::Vector{S}`: (`hist_periods` x 1) vector returning elapsed runtime per period t

"""
function tpf{S<:AbstractFloat}(m::AbstractModel, data::Array{S}, Φ::Function,
                               Ψ::Function, F_ε::Distribution, F_u::Distribution, s_init::Matrix{S};
                               verbose::Symbol=:low, include_presample::Bool=true)
    #--------------------------------------------------------------
    # Set Parameters of Algorithm
    #--------------------------------------------------------------

    # Get tuning parameters from the model
    r_star       = get_setting(m, :tpf_r_star)
    c            = get_setting(m, :tpf_c_star)
    accept_rate  = get_setting(m, :tpf_accept_rate)
    target       = get_setting(m, :tpf_target)
    N_MH         = get_setting(m, :tpf_n_mh_simulations)
    n_particles  = get_setting(m, :tpf_n_particles)
    adaptive     = get_setting(m, :tpf_adaptive)
    xtol         = get_setting(m, :tpf_x_tolerance)
    parallel     = get_setting(m, :use_parallel_workers)

    # Set number of presampling periods
    n_presample_periods = (include_presample) ? 0 : get_setting(m, :n_presample_periods)

    # Initialization of constants and output vectors
    n_observables = size(data,1)
    n_states      = size(s_init)[1]
    T             = size(data,2)
    lik           = zeros(T)
    Neff          = zeros(T)
    times         = zeros(T)

    # Ensuring Φ, Ψ broadcast to matrices
    Φ_bcast(s_t1::Matrix{S}, ε_t1::Matrix{S}) = hcat([Φ(s_t1[:,i], ε_t1[:,i]) for i in 1:size(s_t1)[2]]...)
    Ψ_bcast(s_t1::Matrix{S}, u_t1::Matrix{S}) = hcat([Ψ(s_t1[:,i], u_t1[:,i]) for i in 1:size(s_t1)[2]]...)

    #--------------------------------------------------------------
    # Main Algorithm: Tempered Particle Filter
    #--------------------------------------------------------------

    # Draw initial particles from the distribution of s₀: N(s₀, P₀)
    s_lag_tempered = s_init

    for t = 1:T

        tic()
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("============================================================")
            @show t
        end

        #--------------------------------------------------------------
        # Initialize Algorithm: First Tempering Step
        #--------------------------------------------------------------
        y_t = data[:,t]

        # Remove rows/columns of series with NaN values
        nonmissing          = !isnan(y_t)
        y_t                 = y_t[nonmissing]
        n_observables_t     = length(y_t)
        Ψ_bcast_t           = (x, ε) -> Ψ_bcast(x, ε)[nonmissing, :]
        HH_t                = F_u.Σ.mat[nonmissing, nonmissing]

        # Draw random shock ε
        ε = rand(F_ε, n_particles)

        # Forecast forward one time step
        s_t_nontempered = Φ_bcast(s_lag_tempered, ε)

        # Error for each particle
        p_error = y_t .- Ψ_bcast_t(s_t_nontempered, zeros(n_observables_t, n_particles))

        # Solve for initial tempering parameter φ_1
        if adaptive
            init_Ineff_func(φ) = solve_inefficiency(φ, 2.0*pi, y_t, p_error, HH_t;
                                                    initialize = true) - r_star
            φ_1 = fzero(init_Ineff_func, 1e-30, 1.0, xtol = xtol)
        else
            φ_1 = 0.25
        end

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            @show φ_1
            println("------------------------------")
        end

        # Correct and resample particles
        loglik, id = correction_selection!(φ_1, 0.0, y_t, p_error, HH_t, n_particles,
                                           initialize=true)

        # Update likelihood
        lik[t] += loglik

        # Update arrays for resampled indices
        s_lag_tempered  = s_lag_tempered[:,id]
        s_t_nontempered = s_t_nontempered[:,id]
        ε               = ε[:,id]

        # Tempering initialization
        φ_old = φ_1
        count = 1

        #--------------------------------------------------------------
        # Main Algorithm
        #--------------------------------------------------------------
        while φ_old < 1

            count += 1

            # Get error for all particles
            p_error = y_t .- Ψ_bcast_t(s_t_nontempered, zeros(n_observables_t, n_particles))

            # Define inefficiency function
            init_ineff_func(φ) = solve_inefficiency(φ, φ_old, y_t, p_error, HH_t) - r_star
            fphi_interval = [init_ineff_func(φ_old) init_ineff_func(1.0)]

            # The below boolean checks that a solution exists within interval
            if prod(sign(fphi_interval)) == -1 && adaptive
                φ_new = fzero(init_ineff_func, φ_old, 1., xtol=xtol)
            elseif prod(sign(fphi_interval)) != -1 && adaptive
                φ_new = 1.
            else # fixed φ
                φ_new = 0.5
            end

            if VERBOSITY[verbose] >= VERBOSITY[:low]
                @show φ_new
            end

            # Correct and resample particles
            loglik, id = correction_selection!(φ_new, φ_old, y_t, p_error, HH_t, n_particles)

            # Update likelihood
            lik[t] += loglik

            # Update arrays for resampled indices
            s_lag_tempered  = s_lag_tempered[:,id]
            s_t_nontempered = s_t_nontempered[:,id]
            ε               = ε[:,id]
            p_error         = p_error[:,id]

            # Update value for c
            c = update_c!(c, accept_rate, target)

            if VERBOSITY[verbose] >= VERBOSITY[:low]
                if VERBOSITY[verbose] >= VERBOSITY[:high]
                    @show c
                end
                println("------------------------------")
            end

            # Mutation Step
            accept_vec = zeros(n_particles)
            print("Mutation ")
            tic()

            if parallel
                print("(in parallel) ")
                out = @sync @parallel (hcat) for i=1:n_particles
                    mutation(Φ, Ψ, F_ε, F_u, φ_new, y_t, s_t_nontempered[:,i],
                             s_lag_tempered[:,i], ε[:,i], c, N_MH)
                end
            else
                print("(not parallel) ")
                out = [mutation(Φ, Ψ, F_ε, F_u, φ_new, y_t, s_t_nontempered[:,i],
                                s_lag_tempered[:,i], ε[:,i], c, N_MH)
                       for i = 1:n_particles]
            end
            times[t] = toc()

            for i = 1:n_particles
                s_t_nontempered[:,i] = out[i][1]
                ε[:,i] = out[i][2]
                accept_vec[i] = out[i][3]
            end

            # Calculate average acceptance rate
            accept_rate = mean(accept_vec)

            # Update φ
            φ_old = φ_new
        end

    if VERBOSITY[verbose] >= VERBOSITY[:high]
        println("Out of main while-loop.")
    end

    @show lik[t]
    s_lag_tempered = s_t_nontempered
    print("Completion of one period ")
    toc()
    end

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("=============================================")
    end

    # Return vector of likelihood indexed by time step and Neff
    return Neff[n_presample_periods + 1:end], lik[n_presample_periods + 1:end], times
end


"""
```
get_chol{S<:Float64}(mat::Array{S,2})
```
Calculate and return the Cholesky of a matrix.

"""
@inline function get_chol{S<:Float64}(mat::Array{S,2})
    return Matrix(chol(nearest_spd(mat)))
end

"""
```
update_c!(c_in::Float64, accept_in::Float64, target_in::Float64)
```
Updates value of c by expression that is function of the target and mean acceptance rates.
Returns the new c, in addition to storing it in the model settings.

"""
@inline function update_c!{S<:Float64}(c_in::S, accept_in::S, target_in::S)
    c_out = c_in*(0.95 + 0.1*exp(20*(accept_in - target_in))/(1 + exp(20*(accept_in - target_in))))
    return c_out
end

"""
```
correction_selection!{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1}, p_error::Array{S,2},
    s_lag_tempered::Array{S,2}, ε::Array{S,2}, HH::Array{S,2}, n_particles::Int;
    initialize::Bool=false)
```
Calculate densities, normalize and reset weights, call multinomial resampling, update state and
error vectors, reset error vectors to 1,and calculate new log likelihood.

### Inputs
- `φ_new::S`: current φ
- `φ_old::S`: φ from last tempering iteration
- `y_t::Array{S,1}`: (`n_observables` x 1) vector of observables at time t
- `p_error::Array{S,1}`: A single particle's error: y_t - Ψ(s_t)
- `HH::Array{S,2}`: measurement error covariance matrix, ∑ᵤ
- `n_particles::Int`: number of particles

### Keyword Arguments
- `initialize::Bool`: Flag indicating whether one is solving for incremental weights during
    the initialization of weights; default is `false`.

### Outputs
- `loglik::S`: incremental log likelihood
- `id::Vector{Int}`: vector of indices corresponding to resampled particles
"""
function correction_selection!{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1},
                                   p_error::Array{S,2}, HH::Array{S,2}, n_particles::Int;
                                   initialize::Bool=false)
    # Initialize vector
    incremental_weights = zeros(n_particles)

    # Calculate initial weights
    for n=1:n_particles
        incremental_weights[n] = incremental_weight(φ_new, φ_old, y_t, p_error[:,n], HH,
                                                    initialize=initialize)
    end

    # Normalize weights
    normalized_weights = incremental_weights ./ mean(incremental_weights)

    # Resampling
    id = multinomial_resampling(normalized_weights)

    # Calculate likelihood
    loglik = log(mean(incremental_weights))

    return loglik, id
end

"""
```
incremental_weight{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1}, p_error::Array{S,1},
    HH::Array{S,2}; initialize::Bool=false)
```
### Inputs
- `φ_new::S`: current φ
- `φ_old::S`: φ value before last
- `y_t::Array{S,1}`: Vector of observables for time t
- `p_error::Array{S,1}`: A single particle's error: y_t - Ψ(s_t)
- `HH::Array{S,2}`: Measurement error covariance matrix

### Keyword Arguments
- `initialize::Bool`: Flag indicating whether one is solving for incremental weights during
    the initialization of weights; default is `false`.

### Output
- Returns the incremental weight of single particle
"""
@inline function incremental_weight{S<:Float64}(φ_new::S, φ_old::S, y_t::Array{S,1},
                                       p_error::Array{S,1}, HH::Array{S,2}; initialize::Bool=false)

    # Initialization step (using 2π instead of φ_old)
    if initialize
        return (φ_new/(2*pi))^(length(y_t)/2) * (det(HH)^(-1/2)) *
            exp(-1/2 * p_error' * φ_new * inv(HH) * p_error)[1]

    # Non-initialization step (tempering and final iteration)
    else
        return (φ_new/φ_old)^(length(y_t)/2) *
            exp(-1/2 * p_error' * (φ_new - φ_old) * inv(HH) * p_error)[1]
    end
end

function initialize_function_system{S<:AbstractFloat}(system::System{S})
    # Unpack system
    RRR    = system[:RRR]
    TTT    = system[:TTT]
    HH     = system[:EE] + system[:MM]*system[:QQ]*system[:MM]'
    DD     = system[:DD]
    ZZ     = system[:ZZ]
    QQ     = system[:QQ]
    EE     = system[:EE]
    MM     = system[:MM]
    sqrtS2 = RRR*get_chol(QQ)'

    @inline Φ(s_t1::Vector{S}, ε_t1::Vector{S}) = TTT*s_t1 + sqrtS2*ε_t1
    @inline Ψ(s_t1::Vector{S}, u_t1::Vector{S}) = ZZ*s_t1 + DD + u_t1

    F_ε = Distributions.MvNormal(zeros(size(QQ)[1]), eye(size(QQ)[1]))
    F_u = Distributions.MvNormal(zeros(size(HH)[1]), HH)

    return Φ, Ψ, F_ε, F_u
end

function initialize_state_draws(s0::Vector{Float64}, F_ε::Distribution, Φ::Function,
                                n_parts::Int; burn::Int = 10000, thin::Int = 5)
    s_init = zeros(length(s0), n_parts)
    s_old = s0
    for i in 1:(burn + thin*n_parts)
        ε = rand(F_ε)
        s_new = Φ(s_old, ε)

        if i > burn && i % thin == 0
            draw_index = convert(Int, (i-burn)/thin)
            s_init[:,draw_index] = s_new
        end

        s_old = s_new
    end
    return s_init
end
