using Test, ModelConstructors, DSGE, Dates, FileIO, Random, JLD2

generate_fulldist_forecast_data = false
if VERSION < v"1.5"
    ver = "111"
else
    ver = "150"
end

Random.seed!(1793)

@testset "Test regime switching" begin
    n_reg_temp = 14

    m = Model1002("ss10", custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true)))
    m <= Setting(:date_forecast_start, Date(2020, 6, 30))
    m <= Setting(:date_conditional_end, Date(2020, 6, 30))
    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2020, 3, 31),
                                                3 => Date(2020, 6, 30)))
    m = setup_regime_switching_inds!(m)

    t_s, r_s, c_s = solve(m, regimes = Int[3])

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
    # Which rule to replace with in which periods
    replace_eqcond = Dict{Int, DSGE.EqcondEntry}()
    for i in 3:(n_reg_temp - 1)
        replace_eqcond[i] = DSGE.EqcondEntry(DSGE.ngdp(), [1., 0.])
    end
    replace_eqcond[n_reg_temp] = DSGE.EqcondEntry(AltPolicy(:historical, eqcond, solve), [1., 0.])

    m <= Setting(:regime_eqcond_info, replace_eqcond)

    m <= Setting(:pgap_value, 12.0 : 0.0)
    m <= Setting(:pgap_type, :ngdp)
    m = setup_regime_switching_inds!(m)

    sys_temp = solve(m, regime_switching = true, regimes = collect(1:n_reg_temp),
                     gensys_regimes = UnitRange{Int}[1:2],
                     gensys2_regimes = UnitRange{Int}[2:n_reg_temp])#, apply_altpolicy=false)
    Γ0s = Vector{Matrix{Float64}}(undef, n_reg_temp - 2)
    Γ1s = Vector{Matrix{Float64}}(undef, n_reg_temp - 2)
    Cs = Vector{Vector{Float64}}(undef,  n_reg_temp - 2)
    Ψs = Vector{Matrix{Float64}}(undef, n_reg_temp - 2)
    Πs = Vector{Matrix{Float64}}(undef, n_reg_temp - 2)
    n_endo = length(collect(keys(m.endogenous_states)))
    for i in (3 - 2):(n_reg_temp - 2)
        global Γ0s[i], Γ1s[i], Cs[i], Ψs[i], Πs[i] = haskey(replace_eqcond, i+2) ? replace_eqcond[i+2].alternative_policy.eqcond(m, i+2) : eqcond(m, i + 2)
    end
    Tcal, Rcal, Ccal = DSGE.gensys2(m, Γ0s, Γ1s, Cs, Ψs, Πs, t_s[1:n_endo, 1:n_endo], r_s[1:n_endo, :], c_s[1:n_endo],
                                         T_switch = 12)
    @test Tcal[1] ≈ sys_temp[1][3][1:n_endo, 1:n_endo]

    # Then when start recrusion with rule, should get same TTTs in every perod
    t, r, c = DSGE.ngdp_solve(m, regime_switching = false, regimes = Int[3])

    get_setting(m, :regime_eqcond_info)[n_reg_temp].alternative_policy = DSGE.ngdp()
    Γ0s[end], Γ1s[end], Cs[end], Ψs[end], Πs[end] = eqcond(m, n_reg_temp)
    Tcal, Rcal, Ccal = DSGE.gensys2(m, Γ0s, Γ1s, Cs, Ψs, Πs, t[1:n_endo, 1:n_endo], r[1:n_endo, :], c[1:n_endo])
    for i in 1:length(Tcal)
        @test Tcal[i] ≈ t[1:n_endo, 1:n_endo]
    end
end


@testset "Test regime switching IRFs match Forward Guidance method IRFs" begin
    horizon = 60
    m = Model1002("ss10")
    m <= Setting(:date_forecast_start, Date(2020, 6, 30))
    m <= Setting(:date_conditional_end, Date(2020, 6, 30))
    m <= Setting(:replace_eqcond, false)
    m <= Setting(:gensys2, false)
    m <= Setting(:regime_switching, false)
    m[:ρ_rm] = 0.0

    system = compute_system(m)
    H = 2
    nstates = size(system[:TTT], 1)
    nshocks = size(system[:RRR], 2)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    shocks = m.exogenous_shocks

    states = zeros(nstates, horizon) #, nshocks)
    obs    = zeros(nobs,    horizon) #, nshocks)
    pseudo = zeros(npseudo, horizon) #, nshocks)

    s_0 = zeros(nstates)
    ### IRFs to FFR peg for H periods
    # Set anticipated shocks to peg FFR from t+1 to t+H
    FFRpeg = -m[:Rstarn].value #-0.25/4 #/400

    PsiR1 = 0 #constant
    # ZZ matrix
    PsiR2 = zeros(nstates)
    PsiR2[m.endogenous_states[:R_t]] = 1

    Rht = system[:RRR][:,vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")])]

    bb = zeros(H+1,1)
    MH = zeros(H+1,H+1)
    for hh = 1:H+1
        bb[hh,1] = (FFRpeg - PsiR1 - PsiR2'*(system[:TTT])^hh*s_0)
        MH[hh,:] = PsiR2'*(system[:TTT])^(hh-1)*Rht
    end
    monshocks = MH\bb
    #% Simulate the model
    etpeg = zeros(nshocks,horizon)
    etpeg[vcat(shocks[:rm_sh], shocks[:rm_shl1]:shocks[Symbol("rm_shl$H")]),1] = monshocks
    states[:, :], obs[:, :], pseudo[:, :] = forecast(system, s_0, etpeg, enforce_zlb = false)

    m = Model1002("ss10")
    m <= Setting(:date_forecast_start, Date(2020, 6, 30))
    m <= Setting(:gensys2, true)
    m <= Setting(:replace_eqcond, true)
    replace_eqcond = Dict{Int, DSGE.EqcondEntry}()
    for i in 3:5
        replace_eqcond[i] = DSGE.EqcondEntry(zero_rate(), [1., 0.])
    end
    m <= Setting(:regime_eqcond_info, replace_eqcond)

    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(1990, 3, 31),
                                                3 => Date(2020, 6, 30),
                                                4 => Date(2020, 9, 30),
                                                5 => Date(2020, 12, 31),
                                                6 => Date(2021, 3, 31)))

    m = setup_regime_switching_inds!(m)

    sys = compute_system(m)
    @test !isapprox(sys[4, :TTT], system[:TTT], atol = 1e-5) # Matrix in regime 4 should not match the matrix w/out regime-switching

    s, o, p = forecast(m, sys, zeros(84), zeros(24, 60),  enforce_zlb = false)

    @test all(isapprox.(s[1:20, :], states[1:20, :], atol = 1e-3))
    @test all(isapprox.(o[1:9, :], obs[1:9, :], atol = 1e-3))
    @test all(isapprox.(p[1:10, :], pseudo[1:10, :], atol = 1e-3))
end

@testset "Temporary alternative policies with non-trivial conditional forecasting" begin
    output_vars = [:forecastobs]

    m = Model1002("ss10", custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap =>
                                                                  Setting(:add_altpolicy_pgap, true)))

    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2020, 3, 31),
                                                3 => Date(2020, 6, 30),
                                                4 => Date(2020, 9, 30),
                                                5 => Date(2020, 12, 31),
                                                6 => Date(2021, 3, 31)))
    m <= Setting(:date_conditional_end, Date(2020, 6, 30))
    m <= Setting(:date_forecast_start, Date(2020, 6, 30))
    m <= Setting(:forecast_horizons, 10)
    setup_regime_switching_inds!(m; cond_type = :full)

    fp = dirname(@__FILE__)
    df = load(joinpath(fp, "../reference/regime_switch_data.jld2"), "regime_switch_df_full")
    date_ind = findfirst(df[!, :date] .== Date(2020, 6, 30))
    df[date_ind, :obs_nominalrate] = .3 / 4.
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:regime_eqcond_info, Dict{Int, DSGE.EqcondEntry}(3 => DSGE.EqcondEntry(DSGE.zero_rate(), [1., 0.])))
    m <= Setting(:pgap_type, :ngdp)
    m <= Setting(:pgap_value, 12.0)
    m <= Setting(:gensys2, true)

    # Check zero rate rule causes an error if specified in a conditional period
    m <= Setting(:regime_eqcond_info, Dict{Int, DSGE.EqcondEntry}(4 => DSGE.EqcondEntry(zero_rate(), [1., 0.]), 5 => DSGE.EqcondEntry(zero_rate(), [1., 0.])))
    fcast = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
                                   df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))
    @test fcast[:forecastobs][m.observables[:obs_nominalrate], 1] > .01
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 2]) < 1e-14
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 3]) < 1e-14
    @test all(abs.(fcast[:forecastobs][m.observables[:obs_nominalrate], 4:end]) .> .01)

    m <= Setting(:date_conditional_end, Date(2020, 6, 30))
    m <= Setting(:date_forecast_start, Date(2020, 3, 31))
    fcast = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
                                   df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))
    @test fcast[:forecastobs][m.observables[:obs_nominalrate], 1] > .01
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 2]) > .01
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 3]) < 1e-14
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 4]) < 1e-14
    @test all(abs.(fcast[:forecastobs][m.observables[:obs_nominalrate], 5:end]) .> .01)

    m <= Setting(:date_conditional_end, Date(2020, 3, 31))
    m <= Setting(:date_forecast_start, Date(2020, 3, 31))
    fcast = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
                                   df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))
    @test fcast[:forecastobs][m.observables[:obs_nominalrate], 1] > 0.
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 2]) < 1e-14
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 3]) < 1e-14
    @test abs(fcast[:forecastobs][m.observables[:obs_nominalrate], 4]) < 1e-14
    @test all(abs.(fcast[:forecastobs][m.observables[:obs_nominalrate], 5:end]) .> 0.)
end

@testset "Calculation of regime indices" begin
    m = AnSchorfheide("ss0")

    rss = Vector{Dict{Int, Date}}(undef, 0)
    fss = Vector{Date}(undef, 0)
    ces = Vector{Date}(undef, 0)
    answers = Dict()
    answers[:full] = []
    answers[:none] = []

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2020, 6, 30), 4 => Date(2020, 9, 30)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 6, 30))
    push!(answers[:none], [1:1, 2:60])
    push!(answers[:full], [1:59])

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2020, 6, 30),
                               4 => Date(2021, 12, 31)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 6, 30))
    push!(answers[:none], [1:6, 7:60])
    push!(answers[:full], [1:5, 6:59])

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2020, 12, 31),
                               4 => Date(2021, 12, 31)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 12, 31))
    push!(answers[:none], [1:2, 3:6, 7:60])
    push!(answers[:full], [1:3, 4:57])

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 3, 31), 3 => Date(2021, 12, 31),
                               4 => Date(2022, 12, 31)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 12, 31))
    push!(answers[:none], [1:6, 7:10, 11:60])
    push!(answers[:full], [1:3, 4:7, 8:57])

    push!(rss, Dict{Int, Date}(1 => date_presample_start(m), 2 => Date(2020, 6, 30), 3 => Date(2021, 12, 31),
                               4 => Date(2022, 12, 31)))
    push!(fss, Date(2020, 6, 30))
    push!(ces, Date(2020, 12, 31))
    push!(answers[:none], [1:6, 7:10, 11:60])
    push!(answers[:full], [1:3, 4:7, 8:57])

    for cond_type in [:full, :none]
        for (i, rs, fs, ce) in zip(1:length(rss), rss, fss, ces)
            m <= Setting(:regime_dates, rs)
            m <= Setting(:date_forecast_start, fs)
            m <= Setting(:date_conditional_end, ce)
            setup_regime_switching_inds!(m)
            fcast_reg_inds = DSGE.get_fcast_regime_inds(m, forecast_horizons(m; cond_type = cond_type), cond_type)
            @test fcast_reg_inds == answers[cond_type][i]
        end
    end
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
    m <= Setting(:forecast_horizons, 12)

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
        check_results = load(joinpath(fp, "../reference/regime_switching_fulldist_forecast_version=" * ver * ".jld2"), "fcast_out")
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
            θ_rs[i, 1:n_pvec] = map(x -> (abs(x) < 1e-14) ? 0. : x, θ[i, :])
            j = n_pvec
            for para in DSGE.get_parameters(m)
                if !isempty(para.regimes)
                    for (k, v) in para.regimes[:value]
                        if k > 1
                            j += 1
                            drawv = θ[i, m.keys[para.key]]
                            θ_rs[i, j] = (k  == 2) ? .9 * drawv : drawv
                            if abs(θ_rs[i, j]) < 1e-14
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
                    @test maximum(abs.(check_results[i][cond_type][var] - fcast_out[i][cond_type][var])) < 5e-5
                end
            end
        end
    end

    if generate_fulldist_forecast_data
        fp = dirname(@__FILE__)
        JLD2.jldopen(joinpath(fp, "../reference/regime_switching_fulldist_forecast_version=" * ver * ".jld2"), true, true, true, IOStream) do file
            write(file, "fcast_out", fcast_out)
        end
    end
end

@testset "Test smoothing with regime switching and gensys2 matches plain Kalman filtering and conditional data" begin
    n_reg_temp = 8

    m = Model1002("ss10", custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true)))
    m <= Setting(:date_forecast_start, Date(2020, 3, 31))
    m <= Setting(:date_conditional_end, Date(2020, 3, 31))
    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2019, 12, 31),
                                                3 => Date(2020, 3, 31)))
    m = setup_regime_switching_inds!(m)

    t_s, r_s, c_s = solve(m, regimes = Int[3])

    regime_dates = Dict{Int, Date}()
    regime_dates[1] = date_presample_start(m)
    for (i, date) in zip(2:n_reg_temp, Date(2019,12,31):Dates.Month(3):Date(3000,3,31))
        regime_dates[i] = date
    end
    m <= Setting(:regime_dates, regime_dates)

    m <= Setting(:regime_switching, true)

    m = setup_regime_switching_inds!(m)

    # Use gensy2 for temporary rule
    m <= Setting(:gensys2, true)
    # Replace eqcond with temp rule
    m <= Setting(:replace_eqcond, true)
    # Which rule to replace with in which periods
    replace_eqcond = Dict{Int, DSGE.EqcondEntry}()
    for i in 3:n_reg_temp-1
        replace_eqcond[i] = DSGE.EqcondEntry(DSGE.zero_rate(), [1., 0.])
    end
    replace_eqcond[n_reg_temp] = DSGE.EqcondEntry(AltPolicy(:historical, eqcond, solve), [1., 0.])

    m <= Setting(:regime_eqcond_info, replace_eqcond)

    m <= Setting(:pgap_value, 12.0 : 0.0)
    m <= Setting(:pgap_type, :ngdp)
    m = setup_regime_switching_inds!(m; cond_type = :full)

    sys = compute_system(m)
    df = load(joinpath(dirname(@__FILE__), "../reference/regime_switch_data.jld2"), "regime_switch_df_none")
    df[end, :obs_hours] = NaN
    df[end, :obs_wages] = NaN
    df[end, :obs_consumption] = NaN
    df[end, :obs_nominalrate] = NaN
    m <= Setting(:forecast_smoother, :durbin_koopman)
    histstates, _, _, _ = smooth(m, df, sys; cond_type = :full, draw_states = false)
    kal = DSGE.filter(m, df, sys; cond_type = :full)
    condobs = sys[3, :ZZ] * histstates[:, end] + sys[3, :DD]
    @test histstates[:, end] ≈ kal[:s_T]

    # Now check the implied observables matches conditional data
    @test condobs[1] ≈ df[end, :obs_gdp]
    @test condobs[4] ≈ df[end, :obs_gdpdeflator]
    @test condobs[5] ≈ df[end, :obs_corepce]
    @test condobs[8] ≈ df[end, :obs_investment]
end

nothing
