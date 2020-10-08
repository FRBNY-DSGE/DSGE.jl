using DSGE, HDF5, ModelConstructors, Test

regen = false


m = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true)))
hist_rule = get_setting(m, :alternative_policy)
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
end
