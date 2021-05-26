"""
```
marginal_data_density(m, data; estimation_method = :smc, calculation_method = :incremental_weights)
```

For calculating the log marginal data density for a given posterior sample.

### Inputs

- `m::Union{AbstractDSGEModel,AbstractVARModel}`
- `data::Matrix{Float64}`

### Keyword Arguments

- `estimation_method::Symbol`: either `:smc` or `:mh`
- `calculation_method::Symbol`: either `:incremental_weights` or `:harmonic_mean`
- `parallel::Bool`
- `smc_estimate_file::String`: Specify estimation cloud file to use if estimation_method is SMC
- `bridge_vec::Vector{String}`: Estimation cloud files for all the bridges up to a scratch estimation.
    Used for incremental weights calculation of MDD. Scratch is last.
- `prior_wts::Vector{Float64}`: Weight on the prior for each bridge.
"""
function marginal_data_density(m::Union{AbstractDSGEModel,AbstractVARModel},
                               data::Matrix{Float64} = Matrix{Float64}(undef, 0, 0);
                               estimation_method::Symbol = :smc,
                               calculation_method::Symbol = :incremental_weights,
                               parallel::Bool = false, bridge_vec::Vector = [],
                               smc_estimate_file::String = rawpath(m, "estimate", "smc_cloud.jld2"),
                               prior_wts::Vector{Float64} = zeros(length(bridge_vec)),
                               end_date = Date(3000,1,1))
    if estimation_method == :mh && calculation_method == :incremental_weights
        throw("Can only calculation MDD with incremental weights if the estimation method is :smc")
    end

    if calculation_method == :incremental_weights
        #=if length(prior_wts) > 0 && prior_wts[end] != 0.0
            @warn "Last element of bridges' prior weight assumed to be 0 b/c that's not a bridge estimation."
            prior_wts[end] = 0.0
        end=#
        if length(prior_wts) > 0
            @assert all(prior_wts .>= 0.0)
            @assert length(prior_wts) == length(bridge_vec)
        end

        log_mdd = 0
        prob_ytilde = zeros(BigFloat, length(bridge_vec))
        log_prob_ytilde = zeros(length(bridge_vec))

        for bridge in length(bridge_vec):-1:1 ## Assuming that all bridges are SMC and so incremental method works
            read_h5 = bridge_vec[bridge][end-1:end] == "h5"
            file        = read_h5 ? h5read(bridge_vec[bridge], "smcparams") : load(bridge_vec[bridge])

            cloud, w, W = file["cloud"], file["w"], file["W"]
            n_parts     = sum(W[:, 1])

            w_W = w[:, 2:end] .* W[:, 1:end-1] ./ n_parts
            log_cmdd = sum(log.(sum(w_W, dims = 1))) # sum over particles, take log, sum over param
            cmdd = exp(BigFloat(log_cmdd)) # Need to use non-log version to do below calculation

            log_mdd += log_cmdd
#=
            @show log_cmdd, cmdd
            if bridge == length(bridge_vec)
                prob_ytilde[bridge] = cmdd
                log_prob_ytilde[bridge] = log_cmdd
                #log_mdd += log_cmdd
            else
                @show log(prob_ytilde[bridge+1]), log(prior_wts[bridge+1] * prior_squared + (1.0-prior_wts[bridge+1]) * prob_ytilde[bridge+1]), prior_squared

                prior_adj = prior_wts[bridge+1] * prior_squared + (1.0-prior_wts[bridge+1]) * prob_ytilde[bridge+1]
                log_ans = log_cmdd + log(prior_adj)

                log_prob_ytilde[bridge] = log_ans
                prob_ytilde[bridge] = exp(log_ans)#cmdd / prob_ytilde[bridge+1] *
                    #(prior_wts[bridge+1] * prior_squared + (1.0-prior_wts[bridge+1]) * prob_ytilde[bridge+1])
                #log_mdd += log_ans#log(prob_ytilde[bridge])
            end
            @show log(prob_ytilde[bridge])
=#        end

        read_h5 = smc_estimate_file[end-1:end] == "h5"
        file        = read_h5 ? h5read(smc_estimate_file, "smcparams") : load(smc_estimate_file)
        cloud, w, W = file["cloud"], file["w"], file["W"]
        n_parts     = sum(W[:, 1])

        w_W = w[:, 2:end] .* W[:, 1:end-1] ./ n_parts
        log_mdd += sum(log.(sum(w_W, dims = 1))) # sum over particles, take log, sum over params
#=
        if length(prior_wts) > 0 && prior_wts[1]  > 0 && prob_ytilde[1] < Inf && prob_ytilde[1] > -Inf
            log_mdd += log(prior_wts[1] * prior_squared + (1.0-prior_wts[1]) * prob_ytilde[1])
        end
=#
        return log_mdd
#end
    elseif calculation_method == :harmonic_mean
        free_para_inds = findall(x -> x.fixed == false, get_parameters(m))

        if estimation_method == :smc
            read_h5 = smc_estimate_file[end-1:end] == "h5"
            cloud = load(smc_estimate_file, "cloud")
            params  = get_vals(cloud)
            logpost = get_logpost(cloud)

            return marginal_data_density(params, logpost, free_para_inds)

        elseif estimation_method == :mh
            if isempty(data)
                throw("For calculating the MDD with a :mh sample, one must provide the data as the
                      second positional argument to this function")
            end
            all_params = h5read(rawpath(m, "estimate", "mhsave.h5"), "mhparams")
            all_params = map(Float64, all_params)
            params = Matrix(thin_mh_draws(m, all_params; jstep = 5)')

            n_para = n_parameters(m)
            n_draws = size(params, 2)

            if parallel
                logpost = @distributed (vcat) for i in 1:n_draws
                    posterior!(m, params[:, i], data)
                end
            else
                logpost = Vector{Float64}(undef, n_draws)
                for i in 1:n_draws
                    logpost[i] = posterior!(m, params[:, i], data)
                end
            end

            return marginal_data_density(params, logpost, free_para_inds)

        else
            throw("Invalid estimation method. Must use either :smc or :mh")
        end
    else
        throw("Invalid MDD calculation method. Must use either :incremental_weights or :harmonic_mean")
    end

end

function tt2string(time_temper::Symbol)
    if time_temper == :new
        return "new"
    elseif time_temper == :old
        return "old"
    elseif time_temper == :whole
        return "whole"
    end
end

function marginal_data_density(params::Matrix{Float64}, logpost::Vector{Float64},
                               free_para_inds::Vector{Int64})
    # From margdensim.m
    n_draws = size(params, 2)
    n_free_para = length(free_para_inds)

    θ_all = params
    θ_bar = mean(θ_all, dims = 2)

    p = 0.1:0.1:0.8
    pcrit = map(x -> chisqinvcdf(n_free_para, x), p)
    densfac = mean(logpost)

    θ_free = θ_all[free_para_inds, :]
    θ_bar_free = vec(mean(θ_free, dims = 2))

    # # Computing the covariance matrix with the built-in function
    # Σ_bar = cov(θ_free, 2)
    # Computing covariance matrix the manual way
    Σ_bar = zeros(n_free_para, n_free_para)
    for i in 1:n_draws
        Σ_bar += θ_free[:, i] * θ_free[:, i]'
    end
    Σ_bar = Σ_bar/n_draws - θ_bar_free * θ_bar_free'

    # Outright invert it
    # Σ_bar_inv = inv(Σ_bar)
    # Or
    F = svd(Σ_bar)
    U, S, V = F.U, F.S, F.V
    bigev = findall(x -> x > 1e-6, S)
    parasigdim = length(bigev)
    parasiglndet = 0
    S = diagm(0 => S)
    for i in 1:n_free_para
        if i > parasigdim
            S[i, i] = 0
        else
            parasiglndet = parasiglndet + log(S[i, i])
            S[i, i] = 1/S[i, i]
        end
    end
    Σ_bar_inv = U * S * U'

    all_invlike = Matrix{Float64}(undef, length(p), 0)

    ####################################
    # TEMPORARY
    ####################################
    all_exp_terms = Matrix{Float64}(undef, length(p), n_draws)
    all_indpara = similar(all_exp_terms)
    all_lnfpara = similar(all_exp_terms)
    ####################################

    for i in 1:n_draws
        θ = θ_free[:, i]
        post = logpost[i]
        res = θ - θ_bar_free

        # If outright inverting Σ_bar
        # lnfpara = -log.(p) - .5 * n_free_para * log(2*pi) - .5 * log(det(Σ_bar)) - .5 * res' * Σ_bar_inv * res
        # If not outright inverting Σ_bar
        lnfpara = -log.(p) .- .5 * n_free_para * log(2*pi) .- .5 * parasiglndet .- .5 * res' * Σ_bar_inv * res
        indpara = (res' * Σ_bar_inv * res) .< pcrit
        invlike = exp.(lnfpara .- post .+ densfac) .* indpara

        ####################################
        # TEMPORARY
        ####################################
        all_exp_terms[:, i] = lnfpara .- post .+ densfac

        all_lnfpara[:, i] = lnfpara
        all_indpara[:, i] = indpara
        ####################################

        !any(isinf, invlike) && (all_invlike = hcat(all_invlike, invlike))
    end

    mean_invlike = mean(all_invlike, dims = 2)
    mean_invlike = Base.filter(x -> isfinite(x), mean_invlike)

    return mean(densfac .- log.(mean_invlike))
end

function marginal_data_density_weighted(params::Matrix{Float64},
                                        logpost::Vector{Float64},
                                        free_para_inds::Vector{Int64}, cloud)
    # From margdensim.m
    n_draws = size(params, 2)
    n_free_para = length(free_para_inds)

    θ_all = params
    θ_bar = weighted_mean(cloud)

    p = 0.1:0.1:0.8
    pcrit = map(x -> chisqinvcdf(n_free_para, x), p)
    densfac = mean(logpost)

    θ_free = θ_all[free_para_inds, :]
    θ_bar_free = θ_bar[free_para_inds] #vec(mean(θ_free, 2))

    # # Computing the covariance matrix with the built-in function
    # Σ_bar = cov(θ_free, 2)
    # Computing covariance matrix the manual way
    #Σ_bar = zeros(n_free_para, n_free_para)
    #for i in 1:n_draws
    #    Σ_bar += θ_free[:, i] * θ_free[:, i]'
    #end
    #Σ_bar = Σ_bar/n_draws - θ_bar_free * θ_bar_free'
    # Computing weighted covariance matrix
    Σ_bar = weighted_cov(cloud)
    Σ_bar = Σ_bar[free_para_inds, free_para_inds]

    # Outright invert it
    # Σ_bar_inv = inv(Σ_bar)
    # Or
    F = svd(Σ_bar)
    U, S, V = F.U, F.S, F.V
    bigev = findall(x -> x > 1e-6, S)
    parasigdim = length(bigev)
    parasiglndet = 0
    S = diagm(0 => S)
    for i in 1:n_free_para
        if i > parasigdim
            S[i, i] = 0
        else
            parasiglndet = parasiglndet + log(S[i, i])
            S[i, i] = 1/S[i, i]
        end
    end
    Σ_bar_inv = U * S * U'

    all_invlike = Matrix{Float64}(undef, length(p), 0)

    ####################################
    # TEMPORARY
    ####################################
    all_exp_terms = Matrix{Float64}(undef, length(p), n_draws)
    all_indpara = similar(all_exp_terms)
    all_lnfpara = similar(all_exp_terms)
    ####################################

    for i in 1:n_draws
        θ = θ_free[:, i]
        post = logpost[i]
        res = θ - θ_bar_free

        # If outright inverting Σ_bar
        # lnfpara = -log.(p) - .5 * n_free_para * log(2*pi) - .5 * log(det(Σ_bar)) - .5 * res' * Σ_bar_inv * res
        # If not outright inverting Σ_bar
        lnfpara = -log.(p) - .5 * n_free_para * log(2*pi) - .5 * parasiglndet - .5 * res' * Σ_bar_inv * res
        indpara = (res' * Σ_bar_inv * res) .< pcrit
        invlike = exp.(lnfpara - post + densfac) .* indpara

        ####################################
        # TEMPORARY
        ####################################
        all_exp_terms[:, i] = lnfpara - post + densfac

        all_lnfpara[:, i] = lnfpara
        all_indpara[:, i] = indpara
        ####################################

        !any(isinf, invlike) && (all_invlike = hcat(all_invlike, invlike))
    end

    mean_invlike = mean(all_invlike, dims = 2)
    mean_invlike = Base.filter(x -> isfinite(x), mean_invlike)

    return mean(densfac - log.(mean_invlike))
end


function marginal_data_density_frontier(m::AbstractDSGEModel,
                                        data::Matrix{Float64} = Matrix{Float64}(undef, 0, 0);
                                        estimation_method::Symbol = :smc,
                                        calculation_method = :incremental_weights,
                                        parallel::Bool = false)

    if estimation_method == :mh && calculation_method == :incremental_weights
        throw("Can only calculation the MDD with incremental weights if estimation method is :smc")
    end

    if calculation_method == :incremental_weights
        #file = load(rawpath(m, "estimate", "smc_cloud.jld2", ["adpt="*string(get_setting(m, :tempering_target))]))
        file = load(rawpath(m, "estimate", "smc_cloud.jld2"))
        cloud, w, W = file["cloud"], file["w"], file["W"]
        w_W = w[:, 2:end] .* W[:, 1:end-1]

        return sum(log.(sum(w_W, dims = 1))) # sum across particles, take log, sum across params

    elseif calculation_method == :harmonic_mean
        free_para_inds = findall(x -> x.fixed == false, m.parameters)

        if estimation_method == :smc
            cloud = load(rawpath(m, "estimate", "smc_cloud.jld2"), "cloud")
            params  = get_vals(cloud)
            logpost = get_logpost(cloud)

            return marginal_data_density(params, logpost, free_para_inds)

        elseif estimation_method == :mh
            if isempty(data)
                throw("For calculating the MDD with a :mh sample, one must provide the data as the
                      second positional argument to this function")
            end
            all_params = h5read(rawpath(m, "estimate", "mhsave.h5"), "mhparams")
            all_params = map(Float64, all_params)
            params = thin_mh_draws(m, all_params; jstep = 5)'

            n_para = n_parameters(m)
            n_draws = size(params, 2)

            if parallel
                logpost = @distributed (vcat) for i in 1:n_draws
                    posterior!(m, params[:, i], data)
                end
            else
                logpost = Vector{Float64}(undef, n_draws)
                for i in 1:n_draws
                    logpost[i] = posterior!(m, params[:, i], data)
                end
            end

            return marginal_data_density(params, logpost, free_para_inds)

        else
            throw("Invalid estimation method. Must use either :smc or :mh")
        end
    else
        throw("Invalid MDD calculation method. Must use either " *
              ":incremental_weights or :harmonic_mean")
    end
end
