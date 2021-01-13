function simulate_data(m::AbstractDSGEModel,
                       para_init::OrderedDict{Symbol, Vector{Float64}};
                       filepath::String = "sim_data.h5",
                       burnin::Int = 1, n_periods::Int = 200,
                       s_0::Vector{Float64} = Vector{Float64}(undef, 0),
                       P_0::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
                       meas_error::Float64 = 0.)

    defined_keys = para_init.keys[[isassigned(para_init.keys, i) for i in 1:length(para_init.keys)]]
    defined_vals = [para_init[key] for key in defined_keys]

    # Ensure that para_init is well-formed
    n_para      = length(defined_keys)
    n_para_sets = length(defined_vals[1])
    @assert all(n_para_sets .== map(length, defined_vals[2:end])) "Must give provide the same number of values for each parameter in para_init"

    for i in 1:n_para_sets
        para_init_single = Dict{Symbol, Float64}(zip(defined_keys, [defined_vals[j][i] for j in 1:n_para]))
        simulate_data(m; para_init = para_init_single, filepath = filepath,
                      burnin = burnin, n_periods = n_periods,
                      s_0 = s_0, P_0 = P_0, meas_error = meas_error)
    end
end

function simulate_data(m::AbstractDSGEModel;
                       para_init::OrderedDict{Symbol, Float64} = OrderedDict{Symbol, Float64}(),
                       filepath::String = "sim_data.h5",
                       burnin::Int = 1, n_periods::Int = 200,
                       s_0::Vector{Float64} = Vector{Float64}(undef, 0),
                       P_0::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
                       meas_error::Float64 = 0.)

    # Set initial parameters
    if !isempty(para_init)
        for (key, value) in para_init
            m[key] = value
        end
    end

    answer = ""
    if isfile(filepath)
        invalid_answer = true
        while invalid_answer
            println("WARNING: $filepath already exists. Do you want to overwrite it with a new dataset? (y/n)")
            answer = readline(stdin)
            if (answer == "y" || answer == "n")
                invalid_answer = false
            else
                println("Invalid answer! Must be `y` or `n`")
            end
        end
    end

    # If you want to over-write, or if there is no data currently existing in the $filepath location
    if answer == "y" || !isfile(filepath)
        data = simulate_observables(m, burnin = burnin, n_periods = n_periods, s_0 = s_0, P_0 = P_0, meas_error = meas_error)

        if answer == "y"
            rm(filepath)
            h5write(filepath, "data", data)
        else
            h5write(filepath, "data", data)
        end
        println("Wrote $filepath")

        return data
    elseif answer == "n"
        println("Aborting data simulation because the user chose not to overwrite the existing dataset.")
    end
end

function simulate_observables(m::AbstractDSGEModel;
                              burnin::Int = 1, n_periods::Int = 200,
                              s_0::Vector{Float64} = Vector{Float64}(undef, 0),
                              P_0::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
                              meas_error::Float64 = 1e-3)
    system = compute_system(m; tvis = haskey(get_settings(m), :tvis_information_set))

    n_states = size(system[:TTT], 1)
    n_obs    = length(system[:DD])

    s, ϵ = simulate_states(system[:TTT], system[:RRR], system[:CCC], system[:QQ],
                           burnin = burnin, n_periods = n_periods,
                           s_0 = s_0, P_0 = P_0)

    meas_error_dist = DegenerateMvNormal(zeros(length(system[:DD])),
                                         meas_error*Matrix{Float64}(I, length(system[:DD]), length(system[:DD])))

    y = Matrix{Float64}(undef, n_obs, n_periods)
    for i in 1:n_periods
        y[:, i] = system[:DD] + system[:ZZ]*s[:, i + burnin] + rand(meas_error_dist)
    end

    return y
end

function simulate_states(TTT::Matrix{Float64},
                         RRR::Matrix{Float64},
                         CCC::Vector{Float64},
                         QQ::Matrix{Float64};
                         burnin::Int = 5000,
                         n_periods::Int = 200,
                         ϵ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
                         s_0::Vector{Float64} = Vector{Float64}(undef, 0),
                         P_0::Matrix{Float64} = Matrix{Float64}(undef, 0, 0))

    n_states = size(TTT, 1)
    n_shocks = size(RRR, 2)

    # Initializing shocks
    if isempty(ϵ)
        # Because we parameterize DegenerateMvNormal distributions
        # with the σ not the var
        ϵ = rand(DegenerateMvNormal(zeros(n_shocks), sqrt.(QQ)), burnin + n_periods)
    end

    # Only need shocks starting from t = 1
    # since the states will be 0 indexed, starting from the initial state s_0
    # but will generate one extra shock for the 0 index period and not use
    # it just for convenience of notation
    @assert size(ϵ) == (n_shocks, burnin + n_periods)

    # Initializing states/state covariances
    # Either specify both initial states and state covariances or neither
    @assert (isempty(s_0) && isempty(P_0)) || (!isempty(s_0) && !isempty(P_0))

    if isempty(s_0) && isempty(P_0)
        s_0, P_0 = init_stationary_states(TTT, RRR, CCC, QQ)
    end
    s = Matrix{Float64}(undef, n_states, burnin + n_periods)
    s[:, 1] = s_0

    for i in 2:(burnin + n_periods)
        s[:, i] = iterate_state(s[:, i-1], TTT, RRR, CCC, ϵ[:, i])
    end

    return s, ϵ
end

# Iterates states from s_t to s_{t+1}
function iterate_state(s::Vector{Float64}, TTT::Matrix{Float64}, RRR::Matrix{Float64},
                       CCC::Vector{Float64}, ϵ::Vector{Float64})
    return CCC + TTT*s + RRR*ϵ
end
