# Test that regime-switch in 2020-Q3 to flexible AIT rule runs properly

using Test, ModelConstructors, DSGE, Dates, FileIO, Random

generate_output_data = false

m30 = Model1002("ss30")
m30 <= Setting(:date_forecast_start, Date(2020, 9, 30))
m30 <= Setting(:date_conditional_end, Date(2020, 9, 30))
df = load_data(m30; check_empty_columns = false)
sys30 = compute_system(m30)
m30 <= Setting(:imperfect_awareness_weights, [.5, .5])
sys30_imperf_cred = compute_system(m30)
m30 <= Setting(:imperfect_awareness_weights, [1., 0.])

# Now check m10 generates same results when using flexible AIT policy as an alternative policy
m10 = Model1002("ss10"; custom_settings = Dict{Symbol, Setting}(:add_altpolicy_pgap => Setting(:add_altpolicy_pgap, true),
                                                                :add_altpolicy_ygap => Setting(:add_altpolicy_ygap, true)))
m10 <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m10), 2 => Date(2020, 9, 30)))
m10 <= Setting(:regime_switching, true)
DSGE.setup_regime_switching_inds!(m10)
m10 <= Setting(:alternative_policy, DSGE.flexible_ait())
m10 <= Setting(:pgap_type, :flexible_ait)
m10 <= Setting(:pgap_value, 0.)
m10 <= Setting(:ygap_type, :flexible_ait)
m10 <= Setting(:ygap_value, 12.)
m10 <= Setting(:ait_Thalf, 10.)
m10 <= Setting(:gdp_Thalf, 10.)
m10 <= Setting(:flexible_ait_ρ_smooth, 0.)
m10 <= Setting(:flexible_ait_φ_π, 6.)
m10 <= Setting(:flexible_ait_φ_y, 6.)
sys10_noalt = compute_system(m10; apply_altpolicy = false)
sys10_alt = compute_system(m10; apply_altpolicy = true)

sys30[1, :TTT] ≈ sys10_noalt[1, :TTT]
!(sys30[2, :TTT] ≈ sys10_noalt[2, :TTT])
sys30[1, :TTT] ≈ sys10_alt[1, :TTT]
sys30[2, :TTT] ≈ sys10_alt[2, :TTT]

function param_test(m30, m10; pgap_value = 0.0, ygap_value = 12.0, ait_Thalf = 10.0, gdp_Thalf = 10.0,
                        flexible_ait_ρ_smooth = 0.0, flexible_ait_φ_π = 6.0, flexible_ait_φ_y = 6.0)

    m10 <= Setting(:pgap_value, pgap_value)
    #m10 <= Setting(:ygap_type, :flexible_ait)
    m10 <= Setting(:ygap_value, ygap_value)
    m10 <= Setting(:ait_Thalf, ait_Thalf)
    m10 <= Setting(:gdp_Thalf, gdp_Thalf)
    m10 <= Setting(:flexible_ait_ρ_smooth, flexible_ait_ρ_smooth)
    m10 <= Setting(:flexible_ait_φ_π, flexible_ait_φ_π)
    m10 <= Setting(:flexible_ait_φ_y, flexible_ait_φ_y)
    sys10_noalt = compute_system(m10; apply_altpolicy = false)
    sys10_alt = compute_system(m10; apply_altpolicy = true)

    # ss30
    m30 <= Setting(:pgap_value, pgap_value)
    #m30 <= Setting(:ygap_type, :flexible_ait)
    m30 <= Setting(:ygap_value, ygap_value)
    m30 <= Setting(:ait_Thalf, ait_Thalf)
    m30 <= Setting(:gdp_Thalf, gdp_Thalf)
    m30 <= Setting(:flexible_ait_ρ_smooth, flexible_ait_ρ_smooth)
    m30 <= Setting(:flexible_ait_φ_π, flexible_ait_φ_π)
    m30 <= Setting(:flexible_ait_φ_y, flexible_ait_φ_y)
    sys30 = compute_system(m30)

    @test sys30[1, :TTT] ≈ sys10_noalt[1, :TTT]
    @test !(sys30[2, :TTT] ≈ sys10_noalt[2, :TTT])
    @test sys30[1, :TTT] ≈ sys10_alt[1, :TTT]
    @test sys30[2, :TTT] ≈ sys10_alt[2, :TTT]

    return nothing
end

test_pgap = [0., 2.]
test_ygap = [4., 12.]
test_ait_Thalf = [10.]
test_smooth = [0., 0.5]
test_φ_π = [6., 8.]
test_φ_y = [6., 8.]

@testset "Flexible AIT Policy Change in 2020-Q3" begin
    @simd for pgap in test_pgap
        for ygap in test_ygap
            for ait in test_ait_Thalf
                for smooth in test_smooth
                    for pi in test_φ_π
                        for y in test_φ_y
                            param_test(m30, m10; pgap_value = pgap, ygap_value = ygap, ait_Thalf = ait, gdp_Thalf = ait,
                                       flexible_ait_ρ_smooth = smooth, flexible_ait_φ_π = pi, flexible_ait_φ_y = y)
                        end
                    end
                end
            end
        end
    end

    @test !(sys30_imperf_cred[2, :TTT] ≈ sys30[2, :TTT])
end
