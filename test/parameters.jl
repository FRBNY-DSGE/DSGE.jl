# Tests related to ModelConstructors and making sure
# its code related to parameters work w/DSGE.jl

# subspecs
function sstest(m::AnSchorfheide)

    # Change all the fields of an unfixed parameter
    m <= parameter(:ι_w, 0.000, (0.0, .9999), (0.0,0.9999), ModelConstructors.Untransformed(), Normal(0.0,1.0), fixed=false,
                   description="ι_w: A new parameter.",
                   tex_label="\\iota_w")


    # Change an unfixed parameter to be fixed
    m <= parameter(:ι_p, 0.000, fixed=true,
                   description= "ι_p: Another new parameter",
                   tex_label="\\iota_p")

    # Change a fixed parameter
    m <= parameter(:δ, 0.02,  fixed=true,
                   description="δ: The capital depreciation rate.", tex_label="\\delta" )

    # Overwrite a fixed parameter with an unfixed parameter
    m <= parameter(:ϵ_p, 0.750, (1e-5, 10.),   (1e-5, 10.),     ModelConstructors.Exponential(),    GammaAlt(0.75, 0.4),
                   fixed=false,  scaling = x -> 1 + x/100,
                   description="ϵ_p: No description available.",
                   tex_label="\\varepsilon_{p}")

    steadystate!(m)
end

# Test update! with regime switching
m = Model1002("ss10")
m <= Setting(:regime_switching, true)
m <= Setting(:n_regimes, 2)
m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                            2 => Date(2012, 3, 31)))

ModelConstructors.set_regime_val!(m[:σ_g], 1, 2.523)
ModelConstructors.set_regime_val!(m[:σ_g], 2, 0.; override_bounds = true)
ModelConstructors.set_regime_val!(m[:σ_b], 1, 0.0292)
ModelConstructors.set_regime_val!(m[:σ_b], 2, 0.0292)
ModelConstructors.set_regime_val!(m[:α], 1, 0.1596)
ModelConstructors.set_regime_val!(m[:α], 2, 0.0292)

draw = vec(rand(m.parameters, 1, regime_switching = true))
ModelConstructors.update(m.parameters, draw)
σ_g_ind = findfirst(x->x.key==:σ_g, m.parameters)
σ_b_ind = findfirst(x->x.key==:σ_b, m.parameters)
α_ind = findfirst(x->x.key==:α, m.parameters)
@testset "Test regime switching with two regimes" begin
    @test draw[α_ind] == m.parameters[α_ind].regimes[:value][1]
    @test draw[σ_g_ind] ==  m.parameters[σ_g_ind].regimes[:value][1]
    @test draw[σ_b_ind] ==  m.parameters[σ_b_ind].regimes[:value][1]
    @test draw[end-2] ==  m.parameters[α_ind].regimes[:value][2]
    @test draw[end-1] ==  m.parameters[σ_g_ind].regimes[:value][2]
    @test draw[end] ==  m.parameters[σ_b_ind].regimes[:value][2]
end

m = Model1002("ss10")
m <= Setting(:regime_switching, true)
m <= Setting(:n_regimes, 2)
m <= Setting(:regime_dates, Dict{Int, Date}(1 => date_presample_start(m),
                                            2 => Date(2012, 3, 31),
                                            3 => Date(2015, 3, 31)))

ModelConstructors.set_regime_val!(m[:σ_g], 1, 2.523)
ModelConstructors.set_regime_val!(m[:σ_g], 2, 0.; override_bounds = true)
ModelConstructors.set_regime_val!(m[:σ_g], 3, 2.523)

draw = vec(rand(m.parameters, 1, regime_switching = true))
ModelConstructors.update(m.parameters, draw)

@testset "Test regime switching with two regimes" begin
    @test draw[σ_g_ind] ==  m.parameters[σ_g_ind].regimes[:value][1]
    @test draw[end-1] ==  m.parameters[σ_g_ind].regimes[:value][2]
    @test draw[end] ==  m.parameters[σ_g_ind].regimes[:value][3]
end

m = AnSchorfheide()
sstest(m)

@testset "Test steady-state parameters" begin
    @test m[:ι_w].value == 0.0
    @test m[:ι_w].valuebounds == (0.0, .9999)
    @test m[:ι_w].transform == ModelConstructors.Untransformed()
    @test m[:ι_w].transform_parameterization == (0.0,0.9999)
    @test isa(m[:ι_w].prior.value, Normal)

    @test m[:ι_p].value == 0.0
    @test m[:ι_p].valuebounds == (0.0, 0.0)
    @test isnull(m[:ι_p].prior)
    @test m[:ι_p].fixed == true

    @test m[:δ].value == 0.02
    @test m[:δ].valuebounds == (0.02, 0.02)
    @test isnull(m[:δ].prior)
    @test m[:δ].fixed == true

    @test m[:ϵ_p].value == 0.750
    @test m[:ϵ_p].transform == ModelConstructors.Exponential()
    @test isa(m[:ϵ_p].prior.value, Gamma)
    @test m[:ϵ_p].fixed==false
end

m = AnSchorfheide()
let lastparam = parameter(:p, 0.0)
    for θ in m.parameters
        isa(θ, Parameter) && (lastparam = θ)
    end
    @testset "Check AnSchorfheide last parameter" begin
        @test isa(lastparam, Parameter)
        @test lastparam.value == 0.20*2.237937
    end
end

# transform_to_real_line and transform_to_model_space, acting on the entire parameter vector. they should be inverses!
pvec = m.parameters
vals = transform_to_real_line(pvec)
transform_to_model_space!(m, vals)
@testset "Check parameter transformations for optimization part 2" begin
    @test pvec == m.parameters
end

# all fixed parameters should be unchanged by both transform_to_real_line and transform_to_model_space
@testset "Check fixed parameters are unchanged by optimization transformations" begin
    for θ in m.parameters
        if θ.fixed
            @test θ.value == transform_to_real_line(θ)
            @test θ.value == transform_to_model_space(θ, θ.value)
        end
    end
end
