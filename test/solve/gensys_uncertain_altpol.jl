using DSGE, HDF5, ModelConstructors, Test, Dates

regen = false


m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true)))
hist_rule = DSGE.taylor_rule()
m <= Setting(:alternative_policies, AltPolicy[hist_rule])
m <= Setting(:alternative_policy, DSGE.ngdp())
m <= Setting(:pgap_type, :ngdp)
m <= Setting(:pgap_value, 12.)
T_ngdp, R_ngdp, C_ngdp = get_setting(m, :alternative_policy).solve(m)
T_hist, R_hist, C_hist = hist_rule.solve(m)
rp = joinpath(dirname(@__FILE__), "../reference/")

@testset "Gensys when alternative policies are not credible" begin
    for (i, prob_vec) in enumerate([[1., 0.], [.5, .5], [.999, .001]])
        T_prs, R_prs, C_prs = DSGE.gensys_uncertain_altpol(m, prob_vec, DSGE.AltPolicy[hist_rule],
                                                           apply_altpolicy = true)
        T_prs, R_prs, C_prs = DSGE.augment_states(m, T_prs, R_prs, C_prs)

        m <= Setting(:imperfect_awareness_weights, prob_vec)
        m <= Setting(:uncertain_altpolicy, true)
        TTT, RRR, CCC = solve(m; apply_altpolicy = true) # check automatic calculation
        m <= Setting(:uncertain_altpolicy, false) # to make DSGE.gensys_uncertain_altpol above work

        @test TTT ≈ T_prs
        @test RRR ≈ R_prs
        @test CCC ≈ C_prs

        if regen && i == 2
            h5open(joinpath(rp, "gensys_uncertain_altpol.h5"), "w") do file
                write(file, "T", T_prs)
                write(file, "R", R_prs)
                write(file, "C", C_prs)
            end
        end

        if i == 1
            @test T_prs ≈ T_ngdp
            @test R_prs ≈ R_ngdp
            @test C_prs ≈ C_ngdp
        elseif i == 2
            @test !(T_prs ≈ T_ngdp)
            @test !(R_prs ≈ R_ngdp)
            @test C_prs ≈ C_ngdp    # these should be all zeros
            @test !(T_prs ≈ T_hist)
            @test !(R_prs ≈ R_hist)
            @test C_prs ≈ C_hist    # these should be all zeros
            refpath = joinpath(rp, "gensys_uncertain_altpol.h5")
            @test T_prs ≈ h5read(refpath, "T")
            @test R_prs ≈ h5read(refpath, "R")
            @test C_prs ≈ h5read(refpath, "C")
        else
            @test maximum(abs.(T_prs - T_ngdp)) < 2e-2
            @test maximum(abs.(R_prs - R_ngdp)) < 2e-2
            @test C_prs ≈ C_ngdp    # these should be all zeros
        end
    end

    m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), 2 => Date(2020, 9, 30), 3 => Date(2020, 12, 31),
                                     4 => Date(2021, 3, 31), 5 => Date(2021, 6, 30)))
    m <= Setting(:date_forecast_start, Date(2020, 9, 30))
    m <= Setting(:regime_switching, true)
    setup_regime_switching_inds!(m)

    m <= Setting(:alternative_policies, AltPolicy[hist_rule])
    m <= Setting(:imperfect_awareness_varying_weights,
                 Dict(i => [.5, .5] for i in 2:5))
    m <= Setting(:replace_eqcond, true)
    m <= Setting(:replace_eqcond_func_dict, Dict(i => DSGE.ngdp_replace_eq_entries for i in 2:5))
    m <= Setting(:uncertain_altpolicy, true)
    @test !haskey(DSGE.get_settings(m), :gensys2)
    TTTs, RRRs, CCCs = solve(m; apply_altpolicy = false, regime_switching = true,
                             gensys_regimes = UnitRange{Int}[1:5], # gensys2_regimes = UnitRange{Int}[1:5],
                             regimes = collect(1:5))

    @test !(TTTs[2] ≈ T_ngdp)
    for i in 3:5
        @test TTTs[i] == TTTs[2]
        @test RRRs[i] == RRRs[2]
        @test CCCs[i] == CCCs[2]
    end

    m <= Setting(:imperfect_awareness_varying_weights,
                 Dict(i => [1., 0.] for i in 2:5))
    TTTs2, RRRs2, CCCs2 = solve(m; apply_altpolicy = false, regime_switching = true,
                                gensys_regimes = UnitRange{Int}[1:5], # gensys2_regimes = UnitRange{Int}[1:5],
                                regimes = collect(1:5))

    @test !(TTTs2[2] ≈ TTTs[2])
    for i in 3:5
        @test TTTs2[i] ≈ T_ngdp
        @test RRRs2[i] ≈ R_ngdp
        @test CCCs2[i] ≈ C_ngdp
    end
end
