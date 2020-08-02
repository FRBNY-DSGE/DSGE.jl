file = joinpath(dirname(@__FILE__), "reference/smc_cloud_fix=true.jld2")
cloud = load(file, "cloud")
split_cloud(file, 2)
rejoined_cloud = join_cloud(file, 2)

@testset "Test split and join clouds" begin
    @test cloud.particles           == rejoined_cloud.particles
    @test SMC.get_vals(cloud)       == SMC.get_vals(rejoined_cloud)
    @test SMC.get_loglh(cloud)      == SMC.get_loglh(rejoined_cloud)
    @test SMC.get_old_loglh(cloud)   == SMC.get_old_loglh(rejoined_cloud)
    @test SMC.get_logpost(cloud)    == SMC.get_logpost(rejoined_cloud)
    @test cloud.ESS                 == rejoined_cloud.ESS
    @test cloud.c                   == rejoined_cloud.c
    @test cloud.stage_index         == rejoined_cloud.stage_index
    @test cloud.total_sampling_time == rejoined_cloud.total_sampling_time
    @test cloud.accept              == rejoined_cloud.accept
    @test cloud.n_Φ                 == rejoined_cloud.n_Φ
    @test cloud.resamples           == rejoined_cloud.resamples
    @test cloud.tempering_schedule  == rejoined_cloud.tempering_schedule
end
