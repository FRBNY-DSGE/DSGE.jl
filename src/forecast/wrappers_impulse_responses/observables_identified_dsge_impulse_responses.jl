import DSGE: impulse_responses

function impulse_responses(m::AbstractDSGEModel, paras::Matrix{S},
                           input_type::Symbol, method::Symbol, n_obs_shock::Int,
                           output_vars::Vector{Symbol} =
                           [:irfstates, :irfobs, :irfpseudo]; parallel::Bool = false,
                           permute_mat::Matrix{S} = Matrix{Float64}(undef,0,0),
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false, test_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           do_rev_transform::Bool = false,
                           verbose::Symbol = :high) where {S<:Real}
    if n_obs_shock <= 0
        error("To use method $method, user must specify the index of" *
              " the target observable with keyword `n_obs_shock`.")
    end

    # Set up computation method
    mapfcn = parallel ? pmap : map
    h = impulse_response_horizons(m)
    fcn = if method == :maximum_business_cycle_variance || method == :maxBC
        function _maxBCirf(model, para)
            DSGE.update!(model, para)
            system = compute_system(model)
            states, obs, pseudo = impulse_responses(system, h, frequency_band,
                                                    n_obs_shock, flip_shocks =
                                                    flip_shocks)
            if do_rev_transform
                for (k,v) in m.observable_mappings
                    irf_trans = DSGE.get_irf_transform(v.rev_transform)
                    obs[m.observables[k],:] = irf_trans(obs[m.observables[k],:])
                end
            end
            return states, obs, pseudo
        end
    elseif method in [:cholesky, :choleskyLR, :cholesky_long_run]
        if isempty(permute_mat)
            @warn "Permutation matrix permute_mat is empty. Defaulting to identity."
            permute_mat = Matrix{S}(I, n_observables(m), n_observables(m))
        end
        function _choleskyirf(model, para)
            DSGE.update!(model, para)
            system = compute_system(model)

            obs_shock = zeros(n_observables(model))
            obs_shock[n_obs_shock] = 1.
            states, obs, pseudo = impulse_responses(system, h, permute_mat,
                                                    obs_shock; restriction = method,
                                                    flip_shocks = flip_shocks)
            if do_rev_transform
                for (k,v) in m.observable_mappings
                    irf_trans = DSGE.get_irf_transform(v.rev_transform)
                    obs[m.observables[k],:] = irf_trans(obs[m.observables[k],:])
                end
            end
            return states, obs, pseudo
        end
    elseif method == :regime_switching_structural_irfs
        error("Method $method has not been completely implemented and tested yet.")
        function _struct_irfs(model, para)
            DSGE.update!(model, para)
            regime_system = compute_system(model; regime_switching = true)
            n_regs = n_regimes(regime_system)
            vec_irfstates = Vector{Array{S}}(undef,n_regs)
            vec_irfobs    = Vector{Array{S}}(undef,n_regs)
            vec_irfpseudo = Vector{Array{S}}(undef,n_regs)
            for i = 1:n_regimes(regime_system)
                vec_irfstates[i], vec_irfobs[i], vec_irfpseudo[i] =
                    impulse_responses(m, System(regime_system, i))
            end
            return vec_irfstates, vec_irfobs, vec_irfpseudo
        end
    else
        error("Method $method is not recognized. See the docstring for available impulse response methods.")
    end

    # Compute IRFs
    paras = mapslices(x -> [vec(x)], paras, dims = 2)
    irf_output = mapfcn(para -> fcn(m, para), paras)

    # Reformat output
    states = map(x -> x[1], irf_output)
    obs    = map(x -> x[2], irf_output)
    pseudo = map(x -> x[3], irf_output)

    if create_meansbands
        # Set up metadata and output from IRFs computation
        mb_output_vars = Vector{Matrix{Matrix{S}}}(undef,0)
        class_vars = sort(map(x -> get_class(x), output_vars))
        if :obs in class_vars
            push!(mb_output_vars, obs)
        end
        if :pseudo in class_vars
            push!(mb_output_vars, pseudo)
        end
        if :states in class_vars
            push!(mb_output_vars, states)
        end

        if method != :regime_switching_structural_irfs
            mb_vec = Vector{MeansBands}(undef,0)
            for (mb_i, output_var, class) in zip(1:length(class_vars), mb_output_vars, class_vars)
                metadata = Dict{Symbol,Any}()
                metadata[:para] = input_type
                metadata[:cond_type] = :none
                metadata[:product] = :irf
                metadata[:class] = class
                metadata[:date_inds] = OrderedDict()

                # Set up for loop over variable names
                means = DataFrame()
                bands = Dict{Symbol,DataFrame}()
                names = if class == :states
                    keys(m.endogenous_states)
                elseif class == :obs
                    keys(m.observables)
                elseif class == :pseudo
                    keys(m.pseudo_observables)
                end
                metadata[:indices] = OrderedDict{Symbol,Int}(name => name_i
                                                             for (name_i, name) in enumerate(names))

                # Means and Bands for each variable in a class
                for (name_i,name) in enumerate(names)
                    # Note: obs is Vector{nobs x nperiods} -> for each observable,
                    # so we want to select the IRF of a specific obs, i.e. map(x -> x[obs_index,:]).
                    # This creates a nperiod x ndraws matrix, which we want to transpose
                    # to get a ndraws x nperiod matrix
                    single_var = convert(Matrix{S}, reduce(hcat, map(x -> x[name_i,:], output_var))')
                    means[!,name] = vec(mean(single_var, dims = 1))
                    bands[name]   = find_density_bands(single_var, density_bands;
                                                       minimize = minimize)
                end
                push!(mb_vec, MeansBands(metadata, means, bands))

                # Save MeansBands
                if !test_meansbands
                    tail = method == :cholesky ? :_cholesky : :_maxBC
                    fp = get_meansbands_output_file(m, input_type, :none,
                                                    Symbol(:irf, class, tail),
                                                    forecast_string = forecast_string)
                    dirpath = dirname(fp)
                    isdir(dirpath) || mkpath(dirpath)
                    JLD2.jldopen(fp, true, true, true, IOStream) do file
                        write(file, "mb", mb_vec[mb_i])
                    end
                    println(verbose, :high, "  " * "wrote " * basename(fp))
                end
            end
        end
        return mb_vec
    else
        # Reshape to be nobs x nperiod x ndraw
        if method == :regime_switching_structural_irfs
            n_regs = length(states[1]) # since states is a Vector of Vector{Array}s,
                                       # where the Vector{Array}s is a vector of
                                       # impulse response matrices.
            states = cat(map(x -> [x[i]' for i = 1:n_regs], states)..., dims = 3)
            obs    = cat(map(x -> [x[i]' for i = 1:n_regs], obs)..., dims = 3)
            pseudo = cat(map(x -> [x[i]' for i = 1:n_regs], pseudo)..., dims = 3)
        else
            states = cat(states..., dims = 3)
            obs    = cat(obs..., dims = 3)
            pseudo = cat(pseudo..., dims = 3)
        end
        return states, obs, pseudo
    end
end

function impulse_responses(m::AbstractDSGEModel, paras::Vector{S},
                           input_type::Symbol, method::Symbol, n_obs_shock::Int,
                           output_vars::Vector{Symbol} =
                           [:irfstates, :irfobs, :irfpseudo]; parallel::Bool = false,
                           permute_mat::Matrix{S} = Matrix{Float64}(undef,0,0),
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           create_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           do_rev_transform::Bool = false,
                           verbose::Symbol = :high) where {S<:Real}
    return impulse_responses(m, reshape(paras, 1, length(paras)),
                             input_type, method, n_obs_shock,
                             output_vars; parallel = parallel,
                             permute_mat = permute_mat,
                             frequency_band = frequency_band,
                             flip_shocks = flip_shocks,
                             density_bands = density_bands,
                             create_meansbands = create_meansbands,
                             minimize = minimize,
                             forecast_string = forecast_string,
                             do_rev_transform = do_rev_transform,
                             verbose = verbose)
end
