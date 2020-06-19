using Test, ModelConstructors, DSGE, Dates, FileIO, Random

generate_fulldist_forecast_data = false

Random.seed!(1793)
@testset "Test regime switching" begin
    n_reg_temp = 14

    m = Model1002("ss10", custom_settings = Dict{Symbol, Setting}(:add_pgap => Setting(:add_pgap, true)))

    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                  2 => Date(2020, 3, 31),
                                                  3 => Date(2020, 6, 30)))
    m = setup_regime_switching_inds!(m)

    t_s, r_s, c_s = solve(m, regimes = 3)

    regime_dates = Dict{Int, Date}()
    regime_dates[1] = date_presample_start(m)
    for (i, date) in zip(2:n_reg_temp, Date(2020,3,31):Dates.Month(3):Date(3000,3,31))
        regime_dates[i] = date
    end
    m <= Setting(:regime_dates, regime_dates)

    m <= Setting(:regime_switching, true)

    m = setup_regime_switching_inds!(m)

    # Use gensy2 for temporary rule
    m <= Setting(:gensys2, true)
    # Replace eqcond with temp rule
    m <= Setting(:replace_eqcond, true)
    # Which rule to replace with inn nwhich periods
    replace_eqcond = Dict{Int, Function}()
    for i in 3:n_reg_temp
       replace_eqcond[i] = ngdp_replace_eq_entries
   end

    m <= Setting(:replace_eqcond_func_dict, replace_eqcond)

    m <= Setting(:pgap_value, 12.0 : 0.0)
    m <= Setting(:pgap_type, :ngdp)


    sys_temp = solve(m, regime_switching = true, regimes = 1:n_reg_temp, hist_regimes = 1:2, fcast_regimes = 3:n_reg_temp)
    Γ0s = Vector{Matrix{Float64}}(undef, n_reg_temp-2)
    Γ1s = Vector{Matrix{Float64}}(undef, n_reg_temp-2)
    Cs = Vector{Vector{Float64}}(undef,  n_reg_temp-2)
    Ψs = Vector{Matrix{Float64}}(undef, n_reg_temp-2)
    Πs = Vector{Matrix{Float64}}(undef, n_reg_temp-2)
    n_endo = length(collect(keys(m.endogenous_states)))
    for i in (3-2):(n_reg_temp-2)
        global Γ0s[i], Γ1s[i], Cs[i], Ψs[i], Πs[i] = eqcond(m, i+2, new_policy = true)
    end
    Tcal,Rcal, Ccal = DSGE.gensys_cplus(m, Γ0s, Γ1s, Cs, Ψs, Πs, t_s[1:n_endo, 1:n_endo], r_s[1:n_endo, :], c_s[1:n_endo])
    @test Tcal[1] ≈ sys_temp[1][3][1:n_endo, 1:n_endo]
    # Then when start recrusion with rule, should get same TTTs in every perod
    t, r, c = ngdp_solve(m, regime_switching = false, regimes =3)

    Tcal,Rcal, Ccal = DSGE.gensys_cplus(m, Γ0s, Γ1s, Cs, Ψs, Πs, t[1:n_endo, 1:n_endo], r[1:n_endo, :], c[1:n_endo])
    for i in 1:length(Tcal)
        @test Tcal[i] ≈ t[1:n_endo, 1:n_endo]
    end
##
end

@testset "Test rule switching forecast" begin
    output_vars = [:forecastobs, :histobs, :histpseudo, :forecastpseud]

    m = Model1002("ss10", custom_settings = Dict{Symbol, Setting}(:add_pgap =>
                                                                  Setting(:add_pgap, true)))

    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                  2 => Date(2020, 3, 31),
                                                  3 => Date(2020, 6, 30)))
    m = setup_regime_switching_inds!(m)

    fp = dirname(@__FILE__)
    df = load(joinpath(fp, "../reference/regime_switch_data.jld2"), "regime_switch_df_none")

    m <= Setting(:replace_eqcond, false)
    m <= Setting(:replace_eqcond_func_dict, Dict{Int, Function}(
        3 => ngdp_replace_eq_entries))
    m <= Setting(:pgap_type, :ngdp)
    m <= Setting(:pgap_value, 12.0)

    m <= Setting(:gensys2, false)

    m <= Setting(:alternative_policy, AltPolicy(:ngdp, ngdp_eqcond, ngdp_solve, forecast_init = DSGE.ngdp_forecast_init))

    fcast_altperm = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
                                           df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))

    # Testing ore basic permanent BGDP
    m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_pgap =>
                                                                    Setting(:add_pgap, true)))
    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2020, 3, 31),
                                                3 => Date(2020, 6, 30)))

    m = setup_regime_switching_inds!(m)

    m <= Setting(:pgap_value, 12.0)
    m <= Setting(:pgap_type, :ngdp)

    m <= Setting(:alternative_policy, AltPolicy(:ngdp, ngdp_eqcond, ngdp_solve, forecast_init = DSGE.ngdp_forecast_init))

    fcast_altperm2 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
                                            df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))

    @test fcast_altperm[:forecastobs] == fcast_altperm2[:forecastobs]
end

@testset "Full-distribution forecast with conditioning and regime switching during forecast horizon" begin
    # Now we check regime switching properly draws shocks
    mtest = AnSchorfheide(;testing = true)
    m = Model1002("ss10")
    m <= Setting(:saveroot, saveroot(mtest))
    m <= Setting(:dataroot, dataroot(mtest))
    m <= Setting(:regime_switching, true)
    m <= Setting(:use_parallel_workers, false)
    fp = dirname(@__FILE__)
    overrides = forecast_input_file_overrides(m)
    overrides[:full] = "$(fp)/../reference/mhsave_vint=181115.h5"
    θ = load_draws(m, :full)[1:3, :]
    m <= Setting(:forecast_jstep, 1)
    m <= Setting(:forecast_block_size, 3)

    rss = Vector{Dict{Int, Date}}(undef, 0)
    fss = Vector{Date}(undef, 0)
    ces = Vector{Date}(undef, 0)
    fcast_out = Dict()

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2020, 6, 30), 4 => Date(2020, 9, 30)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 6, 30))

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2020, 6, 30)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 6, 30))

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2020, 12, 31)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 12, 31))

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2021, 12, 31)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 12, 31))

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 6, 30), 3 => Date(2021, 12, 31)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 12, 31))

    output_vars = [:forecastobs, :forecast4qobs, :forecastpseudo, :histpseudo]
    dfs = Dict()
    fp = dirname(@__FILE__)
    dfs[:full] = load(joinpath(fp, "../reference/regime_switch_data.jld2"), "regime_switch_df_full")
    dfs[:semi] = load(joinpath(fp, "../reference/regime_switch_data.jld2"), "regime_switch_df_semi")
    dfs[:none] = load(joinpath(fp, "../reference/regime_switch_data.jld2"), "regime_switch_df_none")

    if !generate_fulldist_forecast_data
        check_results = load(joinpath(fp, "../reference/regime_switching_fulldist_forecast.jld2"), "fcast_out")
    end

    for (i, rs, fs, ce) in zip(1:length(rss), rss, fss, ces)
        # Set up regime switching
        m <= Setting(:regime_dates, rs)
        m <= Setting(:date_forecast_start, fs)
        m <= Setting(:date_conditional_end, ce)
        for i in 1:length(rs)
            adj = (i == 2) ? .9 : 1.
            ModelConstructors.set_regime_val!(m[:σ_g], i, adj * m[:σ_g].value)
            ModelConstructors.set_regime_val!(m[:σ_b], i, adj * m[:σ_b].value)
            ModelConstructors.set_regime_val!(m[:σ_μ], i, adj * m[:σ_μ].value)
            ModelConstructors.set_regime_val!(m[:σ_ztil], i, adj * m[:σ_ztil].value)
            ModelConstructors.set_regime_val!(m[:σ_λ_f], i, adj * m[:σ_λ_f].value)
            ModelConstructors.set_regime_val!(m[:σ_λ_w], i, adj * m[:σ_λ_w].value)
            ModelConstructors.set_regime_val!(m[:σ_r_m], i, adj * m[:σ_r_m].value)
            ModelConstructors.set_regime_val!(m[:σ_σ_ω], i, adj * m[:σ_σ_ω].value)
            ModelConstructors.set_regime_val!(m[:σ_π_star], i, adj * m[:σ_π_star].value)
            ModelConstructors.set_regime_val!(m[:σ_z_p], i, adj * m[:σ_z_p].value)
        end
        setup_regime_switching_inds!(m)

        # Construct parameter matrix
        θ_rs = zeros(size(θ, 1), ModelConstructors.n_parameters_regime_switching(m))
        n_pvec = DSGE.n_parameters(m)
        for i in 1:size(θ_rs, 1)
            θ_rs[i, 1:n_pvec] = map(x -> (abs(x) < 1e-15) ? 0. : x, θ[i, :])
            j = n_pvec
            for para in DSGE.get_parameters(m)
                if !isempty(para.regimes)
                    for (k, v) in para.regimes[:value]
                        if k > 1
                            j += 1
                            drawv = θ[i, m.keys[para.key]]
                            θ_rs[i, j] = (k  == 2) ? .9 * drawv : drawv
                            if abs(θ_rs[i, j]) < 1e-15
                                θ_rs[i, j] = 0.
                            end
                        end
                    end
                end
            end
        end

        fcast_out[i] = Dict()
        if date_conditional_end(m) == Date(2020, 12, 31) && i == 3
            lastrow = deepcopy(dfs[:full][end, :])
            push!(dfs[:full], vcat(Date(2020, 9, 30), lastrow[2:end]...))
            push!(dfs[:full], vcat(Date(2020, 12, 31), lastrow[2:end]...))
            lastrow = deepcopy(dfs[:semi][end, :])
            push!(dfs[:semi], vcat(Date(2020, 9, 30), lastrow[2:end]...))
            push!(dfs[:semi], vcat(Date(2020, 12, 31), lastrow[2:end]...))
        end
        for cond_type in [:full, :none, :semi]
            if i > 2 && cond_type == :none
                continue
            end
            forecast_one(m, :full, cond_type, output_vars, check_empty_columns = false,
                         params = θ_rs, df = dfs[cond_type], verbose = :none)
            output_files = get_forecast_output_files(m, :full, cond_type, output_vars)

            fcast_out[i][cond_type] = Dict()
            for var in keys(output_files)
                fcast_out[i][cond_type][var] = load(output_files[var], "arr")
            end

            if !generate_fulldist_forecast_data
                for var in keys(output_files)
                    @test @test_matrix_approx_eq check_results[i][cond_type][var] fcast_out[i][cond_type][var]
                end
            end
        end
    end

    if generate_fulldist_forecast_data
        fp = dirname(@__FILE__)
        JLD2.jldopen(joinpath(fp, "../reference/regime_switching_fulldist_forecast.jld2"), true, true, true, IOStream) do file
            write(file, "fcast_out", fcast_out)
        end
    end
end

nothing
