path = dirname(@__FILE__())

# Set up arguments
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))

forecast_args = load("$path/../reference/forecast_args.jld2")
df = forecast_args["df"]
system = forecast_args["system"]

# Read expected output
smooth_out = load("$path/../reference/smooth_out.jld2")
exp_states = smooth_out["exp_states"]
exp_shocks = smooth_out["exp_shocks"]
exp_pseudo = smooth_out["exp_pseudo"]

# Smooth without drawing states
states = Dict{Symbol, Matrix{Float64}}()
shocks = Dict{Symbol, Matrix{Float64}}()
pseudo = Dict{Symbol, Matrix{Float64}}()

@testset "Test smoother without drawing states" begin
    for smoother in [:hamilton, :koopman, :carter_kohn, :durbin_koopman]
        m <= Setting(:forecast_smoother, smoother)

        states[smoother], shocks[smoother], pseudo[smoother] =
            smooth(m, df, system; draw_states = false)

        @test @test_matrix_approx_eq exp_states states[smoother]
        @test @test_matrix_approx_eq exp_shocks shocks[smoother]
        @test @test_matrix_approx_eq exp_pseudo pseudo[smoother]
    end
end

# Smooth, drawing states
for smoother in [:carter_kohn, :durbin_koopman]
    m <= Setting(:forecast_smoother, smoother)
    smooth(m, df, system; draw_states = true)
end


# Now check for regime switching
custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_presample_start     => Setting(:date_presample_start, Date(1959, 9, 30)),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_conditional_end, DSGE.quartertodate("2016-Q3")),
    :n_mon_anticipated_shocks => Setting(:n_mon_anticipated_shocks, 6))

m = Model1002("ss10", custom_settings = custom_settings, testing = true)
m <= Setting(:rate_expectations_source, :ois)
# df = load_data(m; check_empty_columns = false, verbose = :none, summary_statistics = :none)
df = load("$path/../reference/regime_switch_data.jld2", "none")

m_rs1 = Model1002("ss10", custom_settings = custom_settings) # pseudo regime switching (no values have second/third regimes)
m_rs1 <= Setting(:rate_expectations_source, :ois)
m_rs1 <= Setting(:regime_switching, true)
m_rs1 <= Setting(:n_regimes, 3)
m_rs1 <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m_rs1), 2 => Date(2010, 3, 31),
                                                3 => Date(2012, 3, 31)))
setup_regime_switching_inds!(m_rs1)

m_rs2 = Model1002("ss51v", custom_settings = custom_settings) # pseudo regime switching (identical values for standard deviations)
m_rs2 <= Setting(:rate_expectations_source, :ois)
m_rs2 <= Setting(:regime_switching, true)
m_rs2 <= Setting(:n_regimes, 3)
m_rs2 <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m_rs1), 2 => Date(2010, 3, 31),
                                                3 => Date(2012, 3, 31)))
setup_regime_switching_inds!(m_rs2)

m_rs3 = Model1002("ss51v", custom_settings = custom_settings) # non-trivial regime switching
m_rs3 <= Setting(:rate_expectations_source, :ois)
m_rs3 <= Setting(:regime_switching, true)
m_rs3 <= Setting(:n_regimes, 3)
m_rs3 <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m_rs1), 2 => Date(2010, 3, 31),
                                                3 => Date(2012, 3, 31)))
setup_regime_switching_inds!(m_rs3)

# Need to set shocks for second and third regimes
prop = 1.
prop3 = .95
for i in 1:3
    if i == 1
        oldprop3 = prop3
        prop3 = prop
    end
    ModelConstructors.set_regime_val!(m_rs2[:σ_g], i, prop * m[:σ_g].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_b], i, prop * m[:σ_b].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_μ], i, prop * m[:σ_μ].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_ztil], i, prop * m[:σ_ztil].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_λ_f], i, prop * m[:σ_λ_f].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_λ_w], i, prop * m[:σ_λ_w].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_r_m], i, prop * m[:σ_r_m].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_σ_ω], i, prop * m[:σ_σ_ω].value)
    ModelConstructors.set_regime_val!(m_rs2[:σ_μ_e], i, prop * m[:σ_μ_e].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_γ], i, prop * m[:σ_γ].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_π_star], i, prop * m[:σ_π_star].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_lr], i, prop * m[:σ_lr].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_z_p], i, prop * m[:σ_z_p].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_tfp], i, prop * m[:σ_tfp].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_gdpdef], i, prop * m[:σ_gdpdef].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_corepce], i, prop * m[:σ_corepce].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_gdp], i, prop * m[:σ_gdp].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs2[:σ_gdi], i, prop * m[:σ_gdi].value, override_bounds = true)

    for j = 1:DSGE.n_mon_anticipated_shocks(m_rs2)
        ModelConstructors.set_regime_val!(m_rs2[Symbol("σ_r_m$(j)")], i, prop * m[Symbol("σ_r_m$(i)")])
    end

    ModelConstructors.set_regime_val!(m_rs3[:σ_g], i, prop3 * m[:σ_g].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_b], i, prop3 * m[:σ_b].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_μ], i, prop3 * m[:σ_μ].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_ztil], i, prop3 * m[:σ_ztil].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_λ_f], i, prop3 * m[:σ_λ_f].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_λ_w], i, prop3 * m[:σ_λ_w].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_r_m], i, prop3 * m[:σ_r_m].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_σ_ω], i, prop3 * m[:σ_σ_ω].value)
    ModelConstructors.set_regime_val!(m_rs3[:σ_μ_e], i, prop3 * m[:σ_μ_e].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_γ], i, prop3 * m[:σ_γ].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_π_star], i, prop3 * m[:σ_π_star].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_lr], i, prop3 * m[:σ_lr].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_z_p], i, prop3 * m[:σ_z_p].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_tfp], i, prop3 * m[:σ_tfp].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_gdpdef], i, prop3 * m[:σ_gdpdef].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_corepce], i, prop3 * m[:σ_corepce].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_gdp], i, prop3 * m[:σ_gdp].value, override_bounds = true)
    ModelConstructors.set_regime_val!(m_rs3[:σ_gdi], i, prop3 * m[:σ_gdi].value, override_bounds = true)

    for j = 1:DSGE.n_mon_anticipated_shocks(m_rs3)
        ModelConstructors.set_regime_val!(m_rs3[Symbol("σ_r_m$(j)")], i, prop3 * m[Symbol("σ_r_m$(i)")])
    end

    if i == 1
        global prop3 = oldprop3
    end
end

states = Dict{Symbol, Matrix{Float64}}()
shocks = Dict{Symbol, Matrix{Float64}}()
pseudo = Dict{Symbol, Matrix{Float64}}()
states_rs1 = Dict{Symbol, Matrix{Float64}}()
shocks_rs1 = Dict{Symbol, Matrix{Float64}}()
pseudo_rs1 = Dict{Symbol, Matrix{Float64}}()
states_rs2 = Dict{Symbol, Matrix{Float64}}()
shocks_rs2 = Dict{Symbol, Matrix{Float64}}()
pseudo_rs2 = Dict{Symbol, Matrix{Float64}}()
states_rs3 = Dict{Symbol, Matrix{Float64}}()
shocks_rs3 = Dict{Symbol, Matrix{Float64}}()
pseudo_rs3 = Dict{Symbol, Matrix{Float64}}()

# Read expected output
exp_states_regime_switch = load("$path/../reference/smooth_out.jld2", "exp_states_regime_switch")
exp_shocks_regime_switch = load("$path/../reference/smooth_out.jld2", "exp_shocks_regime_switch")
exp_pseudo_regime_switch = load("$path/../reference/smooth_out.jld2", "exp_pseudo_regime_switch")
exp_states_regime_switch_draw = load("$path/../reference/smooth_out.jld2", "exp_states_regime_switch_draw") # these are specifically for
exp_shocks_regime_switch_draw = load("$path/../reference/smooth_out.jld2", "exp_shocks_regime_switch_draw") # the Durbin-Koopman smoother
exp_pseudo_regime_switch_draw = load("$path/../reference/smooth_out.jld2", "exp_pseudo_regime_switch_draw")

@testset "Smoothing with regime switching" begin
    system = compute_system(m)
    system_rs1 = compute_system(m_rs1)
    system_rs2 = compute_system(m_rs2)
    system_rs3 = compute_system(m_rs3)
    for smoother in [:hamilton, :koopman, :carter_kohn, :durbin_koopman]
        m <= Setting(:forecast_smoother, smoother)
        m_rs1 <= Setting(:forecast_smoother, smoother)
        m_rs2 <= Setting(:forecast_smoother, smoother)
        m_rs3 <= Setting(:forecast_smoother, smoother)
        states[smoother], shocks[smoother], pseudo[smoother] =
            smooth(m, df, system; draw_states = false)
        states_rs1[smoother], shocks_rs1[smoother], pseudo_rs1[smoother] =
            smooth(m_rs1, df, system_rs1; draw_states = false)
        states_rs2[smoother], shocks_rs2[smoother], pseudo_rs2[smoother] =
            smooth(m_rs2, df, system_rs2; draw_states = false)
        states_rs3[smoother], shocks_rs3[smoother], pseudo_rs3[smoother] =
            smooth(m_rs3, df, system_rs3; draw_states = false)

        @test @test_matrix_approx_eq states_rs1[smoother] states[smoother]
        @test @test_matrix_approx_eq shocks_rs1[smoother] shocks[smoother]
        @test @test_matrix_approx_eq pseudo_rs1[smoother] pseudo[smoother]
        @test @test_matrix_approx_eq states_rs1[smoother] states_rs2[smoother]
        @test @test_matrix_approx_eq shocks_rs1[smoother] shocks_rs2[smoother]
        @test @test_matrix_approx_eq pseudo_rs1[smoother] pseudo_rs2[smoother]
        @test !(states_rs3[smoother] ≈ states_rs2[smoother])
        @test !(shocks_rs3[smoother] ≈ shocks_rs2[smoother])
        @test !(pseudo_rs3[smoother] ≈ pseudo_rs2[smoother])
        if smoother == :durbin_koopman
            @test @test_matrix_approx_eq states_rs3[smoother] exp_states_regime_switch
            @test @test_matrix_approx_eq shocks_rs3[smoother] exp_shocks_regime_switch
            @test @test_matrix_approx_eq pseudo_rs3[smoother] exp_pseudo_regime_switch[1:21, :] # extra rows b/c used more pseudo-obs when generating the test file
        else
            @test maximum(abs.(states_rs3[smoother] - exp_states_regime_switch))  < 5e-1
            @test maximum(abs.(shocks_rs3[smoother] - exp_shocks_regime_switch))  < 2e-1
            @test maximum(abs.(pseudo_rs3[smoother] - exp_pseudo_regime_switch[1:21, :])) < 3e-1
        end
    end

    @info "The following warnings about calling smoothers with draw_states = true are expected."
    for smoother in [:hamilton, :koopman, :carter_kohn, :durbin_koopman]
        m <= Setting(:forecast_smoother, smoother)
        m_rs1 <= Setting(:forecast_smoother, smoother)
        m_rs2 <= Setting(:forecast_smoother, smoother)
        m_rs3 <= Setting(:forecast_smoother, smoother)

        Random.seed!(1793)
        states[smoother], shocks[smoother], pseudo[smoother] =
            smooth(m, df, system; draw_states = true)
        Random.seed!(1793)
        states_rs1[smoother], shocks_rs1[smoother], pseudo_rs1[smoother] =
            smooth(m_rs1, df, system_rs1; draw_states = true)
        Random.seed!(1793)
        states_rs2[smoother], shocks_rs2[smoother], pseudo_rs2[smoother] =
            smooth(m_rs2, df, system_rs2; draw_states = true)
        Random.seed!(1793)
        states_rs3[smoother], shocks_rs3[smoother], pseudo_rs3[smoother] =
            smooth(m_rs3, df, system_rs3; draw_states = true)

        @test @test_matrix_approx_eq states_rs1[smoother] states[smoother]
        @test @test_matrix_approx_eq shocks_rs1[smoother] shocks[smoother]
        @test @test_matrix_approx_eq pseudo_rs1[smoother] pseudo[smoother]
        @test @test_matrix_approx_eq states_rs1[smoother] states_rs2[smoother]
        @test @test_matrix_approx_eq shocks_rs1[smoother] shocks_rs2[smoother]
        @test @test_matrix_approx_eq pseudo_rs1[smoother] pseudo_rs2[smoother]
        @test !(states_rs3[smoother] ≈ states_rs2[smoother])
        @test !(shocks_rs3[smoother] ≈ shocks_rs2[smoother])
        @test !(pseudo_rs3[smoother] ≈ pseudo_rs2[smoother])

        if smoother == :durbin_koopman
            @test maximum(abs.(states_rs3[smoother] - exp_states_regime_switch_draw)) < 3e-5
            # @test @test_matrix_approx_eq states_rs3[smoother] exp_states_regime_switch_draw
            @test @test_matrix_approx_eq shocks_rs3[smoother] exp_shocks_regime_switch_draw
        #    @test @test_matrix_approx_eq pseudo_rs3[smoother] exp_pseudo_regime_switch_draw
        end
    end
end

#=        jldopen("$path/../reference/smooth_out.jld2", "w") do file
            file["exp_pseudo_regime_switch"] = pseudo_rs3[:durbin_koopman]
            file["exp_states_regime_switch"] = states_rs3[:durbin_koopman]
            file["exp_shocks_regime_switch"] = shocks_rs3[:durbin_koopman]
            file["exp_states"] = exp_states
            file["exp_shocks"] = exp_shocks
            file["exp_pseudo"] = exp_pseudo
            file["exp_pseudo_regime_switch_draw"] = exp_pseudo_regime_switch_draw
            file["exp_states_regime_switch_draw"] = exp_states_regime_switch_draw
            file["exp_shocks_regime_switch_draw"] = exp_shocks_regime_switch_draw
        end =#

nothing
