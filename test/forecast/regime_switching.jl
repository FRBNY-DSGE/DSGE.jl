using Test, ModelConstructors, DSGE, Dates

@testset "Test regime switching" begin
    n_reg_temp = 14

    m = Model1002("ss60", custom_settings = Dict{Symbol, Setting}(:add_pgap => Setting(:add_pgap, true)))

    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                  2 => Date(2020, 3, 31),
                                                  3 => Date(2020, 6, 30)))
    m = setup_regime_switching_inds(m)

    t_s, r_s, c_s = solve(m, regimes=3)

    regime_dates = Dict{Int, Date}()
    regime_dates[1] = date_presample_start(m)
    for (i, date) in zip(2:n_reg_temp, Date(2020,3,31):Dates.Month(3):Date(3000,3,31))
        regime_dates[i] = date
    end
    m <= Setting(:regime_dates, regime_dates)

    m <= Setting(:regime_switching, true)

    m = setup_regime_switching_inds(m)

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

    m = Model1002("ss60", custom_settings = Dict{Symbol, Setting}(:add_pgap =>
                                                                  Setting(:add_pgap, true)))
    df = load_data(m, try_disk =false, check_empty_columns = false)

    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                  2 => Date(2020, 3, 31),
                                                  3 => Date(2020, 6, 30)))
    m = setup_regime_switching_inds(m)

    m <= Setting(:replace_eqcond, false)
    m <= Setting(:replace_eqcond_func_dict, Dict{Int, Function}(
        3 => ngdp_replace_eq_entries))
    m <= Setting(:pgap_type, :ngdp)
    m <= Setting(:pgap_value, 12.0)

    m <= Setting(:gensys2, false)

    m <= Setting(:alternative_policy, AltPolicy(:ngdp, ngdp_eqcond, ngdp_solve, forecast_init = DSGE.ngdp_forecast_init))

    fcast_altperm = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
                                           df; regime_switching = true, n_regimes = get_setting(m, :n_regimes),
                                           )
    # Testing ore basic permanent BGDP
    m = Model1002("ss60"; custom_settings = Dict{Symbol, Setting}(:add_pgap =>
                                                                    Setting(:add_pgap, true)))
    m <= Setting(:regime_switching, true)
    m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                                2 => Date(2020, 3, 31),
                                                3 => Date(2020, 6, 30)))

    m = setup_regime_switching_inds(m)

    m <= Setting(:pgap_value, 12.0)
    m <= Setting(:pgap_type, :ngdp)

    m <= Setting(:alternative_policy, AltPolicy(:ngdp, ngdp_eqcond, ngdp_solve, forecast_init = DSGE.ngdp_forecast_init))

    fcast_altperm2 = DSGE.forecast_one_draw(m, :mode, :full, output_vars, map(x -> x.value, m.parameters),
                                            df; regime_switching = true, n_regimes = get_setting(m, :n_regimes))

    @test fcast_altperm[:forecastobs] == fcast_altperm2[:forecastobs]
end
