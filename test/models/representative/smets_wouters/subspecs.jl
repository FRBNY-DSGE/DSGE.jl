@testset "subspecs: ss0 default subspec" begin
    m = SmetsWouters("ss0")
    @test subspec(m) == "ss0"
end

@testset "subspecs: ss1" begin
    m = SmetsWouters("ss1")
    @test subspec(m) == "ss1"
end

@testset "subspecs: ss2 widens beta parameter bounds to (0.0, 1.0)" begin
    m = SmetsWouters("ss2")
    @test subspec(m) == "ss2"

    for key in (:ζ_p, :ι_p, :h, :ppsi, :ζ_w, :ι_w, :ρ, :ρ_g, :ρ_b, :ρ_μ, :ρ_z, :ρ_λ_f, :ρ_λ_w, :ρ_rm, :η_gz, :η_λ_f, :η_λ_w)
        @test m[key].valuebounds == (0.0, 1.0)
    end

    @test m[:ζ_p].value  ≈ 0.7813
    @test m[:h].value    ≈ 0.7205
    @test m[:ρ_g].value  ≈ 0.9930
    @test m[:ρ_rm].value ≈ 0.3000
end

@testset "subspecs: ss3 builds on ss2 and fixes ρ_g" begin
    m = SmetsWouters("ss3")
    @test subspec(m) == "ss3"

    @test m[:ζ_p].valuebounds == (0.0, 1.0)

    @test m[:ρ_g].fixed == true
    @test m[:ρ_g].value ≈ 0.9930
end

@testset "subspecs: ss4 sets measurement error parameters" begin
    m = SmetsWouters("ss4")
    @test subspec(m) == "ss4"

    scale = 1/5
    @test m[:e_y].value ≈ scale * 0.868241996
    @test m[:e_L].value ≈ scale * 0.283703727
    @test m[:e_w].value ≈ scale * 0.600015046
    @test m[:e_π].value ≈ scale * 0.600386299
    @test m[:e_R].value ≈ scale * 0.846722459
    @test m[:e_c].value ≈ scale * 0.710297826
    @test m[:e_i].value ≈ scale * 2.405946712

    for key in (:e_y, :e_L, :e_w, :e_π, :e_R, :e_c, :e_i)
        @test m[key].fixed == true
    end
end

@testset "subspecs: ss5 sets Lmean and fixes γ" begin
    m = SmetsWouters("ss5")
    @test subspec(m) == "ss5"

    @test m[:Lmean].value ≈ 0.0
    @test m[:Lmean].fixed == false

    @test m[:γ].fixed == true
    @test m[:γ].value ≈ 0.4312
end

@testset "subspecs: ss6 uses diffuse priors" begin
    m = SmetsWouters("ss6")
    @test subspec(m) == "ss6"

    @test m[:α].value      ≈ 0.24
    @test m[:α].fixed      == false
    @test m[:δ].fixed      == true
    @test m[:Upsilon].fixed == true
    @test m[:λ_w].fixed    == true
    @test m[:ϵ_p].fixed    == true
    @test m[:ϵ_w].fixed    == true
    @test m[:g_star].fixed  == true
    @test m[:γ].fixed      == false
    @test m[:γ].value      ≈ 0.3982
    @test m[:Lmean].value  ≈ 875.0
end

@testset "subspecs: invalid subspec throws error" begin
    @test_throws ErrorException SmetsWouters("ss99")
end

nothing
