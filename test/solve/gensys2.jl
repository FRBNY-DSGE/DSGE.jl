using DSGE, ModelConstructors, Dates, Test

@testset "sys1: rule-switch (1st forecast quarter), no parameter-switch. sys2: rule-switch (1st forecast quarter), no parameter-switch (parameter switch with same values)" begin
    m = Model1002("ss10")

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:α], 1, m[:α].value) #2.523)
    ModelConstructors.set_regime_val!(m[:α], 2, m[:α].value)
    ModelConstructors.set_regime_val!(m[:α], 3, m[:α].value) #
    ModelConstructors.set_regime_val!(m[:α], 4, m[:α].value)
    sys2 = compute_system(m)

    @test sys1[1, :TTT] == sys2[1, :TTT]
    @test sys1[2, :TTT] == sys2[2, :TTT]
    @test sys1[3, :TTT] == sys2[3, :TTT]
    @test sys1[4, :TTT] == sys2[4, :TTT]
end


@testset "sys1: no rule-swtich, sys2 rule-switch in 1st forecast quarter (so system 3 should be different)" begin
    m = Model1002("ss10")
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter
    ModelConstructors.set_regime_val!(m[:α], 1, m[:α].value) #2.523)
    ModelConstructors.set_regime_val!(m[:α], 2, m[:α].value)
    ModelConstructors.set_regime_val!(m[:α], 3, m[:α].value) #
    ModelConstructors.set_regime_val!(m[:α], 4, m[:α].value)

    sys1 = compute_system(m)

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)

    sys2 = compute_system(m)

    @test sys1[1, :TTT] == sys2[1, :TTT]
    @test sys1[2, :TTT] == sys2[2, :TTT]
    @test sys1[3, :TTT] != sys2[3, :TTT]
    @test sys1[4, :TTT] == sys2[4, :TTT]
end


@testset "sys1: no rule-switch, no parameter-switch. sys2: rule- and parameter-switch in 1st forecast quarter (so system 3 should be different)" begin
    m = Model1002("ss10")
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:α], 1, m[:α].value) #2.523)
    ModelConstructors.set_regime_val!(m[:α], 2, m[:α].value)
    ModelConstructors.set_regime_val!(m[:α], 3, 0.00001) #
    ModelConstructors.set_regime_val!(m[:α], 4, m[:α].value)

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)

    sys2 = compute_system(m)

    @test sys1[1, :TTT] == sys2[1, :TTT]
    @test sys1[2, :TTT] == sys2[2, :TTT]
    @test sys1[3, :TTT] != sys2[3, :TTT]
    @test sys1[4, :TTT] == sys2[4, :TTT]
end

@testset "sys1: no rule-switch, no parameter-switch. sys2: rule-switch in 1st forecast quarter, parameter-switch in 2nd forecastquarter" begin
    m = Model1002("ss10")
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:α], 1, m[:α].value) #2.523)
    ModelConstructors.set_regime_val!(m[:α], 2, m[:α].value)
    ModelConstructors.set_regime_val!(m[:α], 3, m[:α].value)
    ModelConstructors.set_regime_val!(m[:α], 4, 0.00001) #

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)

    sys2 = compute_system(m)

    @test sys1[1, :TTT] == sys2[1, :TTT]
    @test sys1[2, :TTT] == sys2[2, :TTT]
    @test sys1[3, :TTT] != sys2[3, :TTT]
    @test sys1[4, :TTT] != sys2[4, :TTT]
end

@testset "sys1: no rule-switch, no parameter-swtich. sys2: no rule-switch, yes parameter-switch (in 1st forecast quarter)" begin
    m = Model1002("ss10")
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:α], 1, m[:α].value) #2.523)
    ModelConstructors.set_regime_val!(m[:α], 2, m[:α].value)
    ModelConstructors.set_regime_val!(m[:α], 3, 0.00001) #
    ModelConstructors.set_regime_val!(m[:α], 4, m[:α].value)

    sys2 = compute_system(m)

    @test sys1[1, :TTT] == sys2[1, :TTT]
    @test sys1[2, :TTT] == sys2[2, :TTT]
    @test sys1[3, :TTT] != sys2[3, :TTT]
    @test sys1[4, :TTT] == sys2[4, :TTT]

end


### Now with QQ
@testset "sys1: rule-switch (1st forecast quarter), no parameter-switch. sys2: rule-switch (1st forecast quarter), no parameter-switch (parameter switch with same values)" begin
    m = Model1002("ss10")

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:σ_g], 1, m[:σ_g].value) #2.523)
    ModelConstructors.set_regime_val!(m[:σ_g], 2, m[:σ_g].value)
    ModelConstructors.set_regime_val!(m[:σ_g], 3, m[:σ_g].value) #
    ModelConstructors.set_regime_val!(m[:σ_g], 4, m[:σ_g].value)
    sys2 = compute_system(m)

    @test sys1[1, :QQ] == sys2[1, :QQ]
    @test sys1[2, :QQ] == sys2[2, :QQ]
    @test sys1[3, :QQ] == sys2[3, :QQ]
    @test sys1[4, :QQ] == sys2[4, :QQ]
end

@testset "sys1: no rule-swtich, sys2 rule-switch in 1st forecast quarter (so system 3 should be different)" begin
    # rule switching doesn't affect QQ here
    m = Model1002("ss10")

    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:σ_g], 1, m[:σ_g].value) #2.523)
    ModelConstructors.set_regime_val!(m[:σ_g], 2, m[:σ_g].value)
    ModelConstructors.set_regime_val!(m[:σ_g], 3, m[:σ_g].value) #
    ModelConstructors.set_regime_val!(m[:σ_g], 4, m[:σ_g].value)

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)

    sys2 = compute_system(m)

    @test sys1[1, :QQ] == sys2[1, :QQ]
    @test sys1[2, :QQ] == sys2[2, :QQ]
    @test sys1[3, :QQ] == sys2[3, :QQ]
    @test sys1[4, :QQ] == sys2[4, :QQ]
end


@testset "sys1: no rule-switch, no parameter-switch. sys2: rule- and parameter-switch in 1st forecast quarter (so system 3 should be different)" begin
    m = Model1002("ss10")

    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:σ_g], 1, m[:σ_g].value) #2.523)
    ModelConstructors.set_regime_val!(m[:σ_g], 2, m[:σ_g].value)
    ModelConstructors.set_regime_val!(m[:σ_g], 3, 0.00001) #
    ModelConstructors.set_regime_val!(m[:σ_g], 4, m[:σ_g].value)

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)

    sys2 = compute_system(m)

    @test sys1[1, :QQ] == sys2[1, :QQ]
    @test sys1[2, :QQ] == sys2[2, :QQ]
    @test sys1[3, :QQ] != sys2[3, :QQ]
    @test sys1[4, :QQ] == sys2[4, :QQ]
end

@testset "sys1: no rule-switch, no parameter-switch. sys2: rule-switch in 1st forecast quarter, parameter-switch in 2nd forecastquarter" begin
    # rule switching doesn't affect QQ here
    m = Model1002("ss10")

    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:σ_g], 1, m[:σ_g].value) #2.523)
    ModelConstructors.set_regime_val!(m[:σ_g], 2, m[:σ_g].value)
    ModelConstructors.set_regime_val!(m[:σ_g], 3, m[:σ_g].value) #
    ModelConstructors.set_regime_val!(m[:σ_g], 4, 0.00001)

    # Use gensys2 for anticipated policy rule change
    m <= Setting(:gensys2, true)
    # Number of periods rule is in effect
    m <= Setting(:n_rule_periods, 1)


    sys2 = compute_system(m)

    @test sys1[1, :QQ] == sys2[1, :QQ]
    @test sys1[2, :QQ] == sys2[2, :QQ]
    @test sys1[3, :QQ] == sys2[3, :QQ]
    @test sys1[4, :QQ] != sys2[4, :QQ]
end


@testset "sys1: no rule-switch, no parameter-swtich. sys2: no rule-switch, yes parameter-switch (in 1st forecast quarter)" begin
    m = Model1002("ss10")
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)
    # Number of regimes OVERALL
    m <= Setting(:n_regimes, 4)
    # Number of regimes in forcast
    m <= Setting(:n_fcast_regimes, 2)
    # Number of regimes in history
    m <= Setting(:n_hist_regimes, 2)
    # Dates of all regime switches (including both rule and parameter switches)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30))) #Parameter and/or rule switch in 2nd forecast quarter

    sys1 = compute_system(m)

    ModelConstructors.set_regime_val!(m[:σ_g], 1, m[:σ_g].value) #2.523)
    ModelConstructors.set_regime_val!(m[:σ_g], 2, m[:σ_g].value)
    ModelConstructors.set_regime_val!(m[:σ_g], 3, 0.00001) #
    ModelConstructors.set_regime_val!(m[:σ_g], 4, m[:σ_g].value)

    sys2 = compute_system(m)

    @test sys1[1, :QQ] == sys2[1, :QQ]
    @test sys1[2, :QQ] == sys2[2, :QQ]
    @test sys1[3, :QQ] != sys2[3, :QQ]
    @test sys1[4, :QQ] == sys2[4, :QQ]
end

@testset "Adding redundant regimes does not affect gensys2" begin
    # Add additional check for using replace_eqcond "falsely"
    m = Model1002("ss10")
    # Always need regime switching on, even if you're not switching parameters
    m <= Setting(:regime_switching, true)

    m <= Setting(:date_forecast_start, Date(2020, 9, 30))
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30),
                                     5 => Date(2020, 12, 31),
                                     6 => Date(2021, 3, 31))) #Parameter and/or rule switch in 2nd forecast quarter
    setup_regime_switching_inds!(m)

    m <= Setting(:gensys2, true)
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:regime_eqcond_info, Dict(4 => DSGE.EqcondEntry(DSGE.zero_rate(), [1., 0.])))
    sys1 = compute_system(m)
    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), #Start
                                     2 => Date(1990, 3, 31), # Parameter switch in history
                                     3 => Date(2020, 6, 30), # Parameter and/or rule switch in 1st forecast quarter
                                     4 => Date(2020, 9, 30),
                                     5 => Date(2020, 12, 31),
                                     6 => Date(2021, 3, 31),
                                     7 => Date(2021, 6, 30))) #Parameter and/or rule switch in 2nd forecast quarter
    setup_regime_switching_inds!(m)
    sys2 = compute_system(m)

    for i in 1:n_regimes(sys1)
        @test sys1[i, :TTT] ≈ sys2[i, :TTT]
        @test sys1[i, :RRR] ≈ sys2[i, :RRR]
        @test sys1[i, :CCC] ≈ sys2[i, :CCC]
    end
    @test sys2[6, :TTT] ≈ sys[7, :TTT]
    @test sys2[6, :RRR] ≈ sys[7, :RRR]
    @test sys2[6, :CCC] ≈ sys[7, :CCC]
end
