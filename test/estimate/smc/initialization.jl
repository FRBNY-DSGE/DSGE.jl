write_test_output = false

path = dirname(@__FILE__)

###################################################################
# Test: initial_draw!()
###################################################################
m = AnSchorfheide()

function loglik_fn(para::ModelConstructors.ParameterVector, data::Matrix{Float64})::Float64
    DSGE.update!(m, para)
    DSGE.likelihood(m, data; sampler = false, catch_errors = true,
               use_chand_recursion = get_setting(m, :use_chand_recursion), verbose = :low)
end

save = normpath(joinpath(dirname(@__FILE__),"save"))
m <= Setting(:saveroot, save)

data = h5read(joinpath(path, "reference/smc.h5"), "data")

m <= Setting(:n_particles, 400)
m <= Setting(:n_Φ, 100)
m <= Setting(:λ, 2.0)
m <= Setting(:n_smc_blocks, 1)
m <= Setting(:use_parallel_workers, false)
m <= Setting(:step_size_smc, 0.5)
m <= Setting(:n_mh_steps_smc, 1)
m <= Setting(:resampler_smc, :polyalgo)
m <= Setting(:target_accept, 0.25)

m <= Setting(:mixture_proportion, .9)
m <= Setting(:tempering_target, 0.95)
m <= Setting(:resampling_threshold, .5)
m <= Setting(:use_fixed_schedule, true)

####################################################################
init_cloud = SMC.Cloud(length(m.parameters), get_setting(m, :n_particles))

@everywhere Random.seed!(42)
SMC.initial_draw!(loglik_fn, m.parameters, data, init_cloud)

if write_test_output
    JLD2.jldopen(joinpath(path, "reference/initial_draw_out.jld2"), "w") do file
        write(file, "cloud", init_cloud)
    end
end

saved_init_cloud = load(joinpath(path, "reference/initial_draw_out.jld2"), "cloud")

@testset "Initial draw" begin
    @test @test_matrix_approx_eq SMC.get_vals(init_cloud) SMC.get_vals(saved_init_cloud)
    @test @test_matrix_approx_eq SMC.get_loglh(init_cloud) SMC.get_loglh(saved_init_cloud)
end

###################################################################
# Test: one_draw()
###################################################################
#=
# Commented out b/c difficulty reconstructing parameters. This feature is
# already tested in SMC.jl, so it should still work anyway.
m = AnSchorfheide()

file = JLD2.jldopen(joinpath(path, "reference/one_draw_in.jld2"), "r")
    data = file["data"]
close(file)
file = JLD2.jldopen(joinpath(path, "reference/one_draw_in_parameter_fields.jld2"), "r")
    parameter_fields = file["parameter_fields"]
close(file)
p = parameter_fields
parameters = map(x -> (:scaling in fieldnames(typeof(p[x]))) ?
                 parameter(p[x][:key], p[x][:value], p[x][:valuebounds],
                           p[x][:transform_parameterization], p[x][:transform],
                           p[x][:prior]; fixed = p[x][:fixed],
                           description = p[x][:description],
                           tex_label = p[x][:tex_label],
                           scaling = p[x][:scaling]) :
                 parameter(p[x][:key], p[x][:value], p[x][:valuebounds],
                           p[x][:transform_parameterization], p[x][:transform],
                           p[x][:prior]; fixed = p[x][:fixed],
                           description = p[x][:description],
                           tex_label = p[x][:tex_label]), 1:length(p))
parameters = convert(Vector{AbstractParameter{Float64}}, parameters)

draw = SMC.one_draw(loglik_fn, parameters, data)

if write_test_output
    JLD2.jldopen(joinpath(path, "reference/one_draw_out.jld2"), true, true, true, IOStream) do file
        file["draw"] = draw
    end
end

test_draw = JLD2.jldopen(joinpath(path, "reference/one_draw_out.jld2"), "r") do file
    file["draw"]
end

###################################################################
@testset "One draw" begin
    @test draw[1] == test_draw[1]
    @test draw[2] ≈ test_draw[2]
end
=#
###################################################################
# Test: draw_likelihood()
###################################################################
draw = JLD2.jldopen(joinpath(path, "reference/one_draw_out.jld2"), "r") do file
    file["draw"]
end

draw_lik = SMC.draw_likelihood(loglik_fn, parameters, data, vec(draw[1]))

if write_test_output
    JLD2.jldopen(joinpath(path, "reference/draw_likelihood_out.jld2"), true, true, true, IOStream) do file
        file["draw_lik"] = draw_lik
    end
end
test_draw_lik = JLD2.jldopen(joinpath(path, "reference/draw_likelihood_out.jld2"), "r") do file
    file["draw_lik"]
end

###################################################################
@testset "Draw Likelihood" begin
    @test draw_lik[1][1] ≈ test_draw_lik[1][1]
    @test draw_lik[2] == test_draw_lik[2]
end

###################################################################
# Test: initialize_likelihoods!()
###################################################################
SMC.initialize_likelihoods!(loglik_fn, parameters, data, init_cloud)

if write_test_output
    JLD2.jldopen(joinpath(path, "reference/initialize_likelihood_out.jld2"), true, true, true, IOStream) do file
        file["init_lik_cloud"] = init_cloud
    end
end
test_init_cloud = JLD2.jldopen(joinpath(path, "reference/initialize_likelihood_out.jld2"), "r") do file
    file["init_lik_cloud"]
end

###################################################################
@testset "Initialize Likelihoods" begin
    @test @test_matrix_approx_eq SMC.get_vals(init_cloud) SMC.get_vals(test_init_cloud)
    @test @test_matrix_approx_eq SMC.get_loglh(init_cloud) SMC.get_loglh(test_init_cloud)
end


###################################################################
# Test: initialize_cloud_settings!()
###################################################################
SMC.initialize_cloud_settings!(init_cloud)

if write_test_output
    JLD2.jldopen(joinpath(path, "reference/initialize_cloud_settings.jld2"), true, true, true, IOStream) do file
        file["init_cloud"] = init_cloud
    end
end
test_init_cloud = JLD2.jldopen(joinpath(path, "reference/initialize_cloud_settings.jld2"), "r") do file
    file["init_cloud"]
end

###################################################################
@testset "Initialize Cloud" begin
    @test init_cloud.ESS         == test_init_cloud.ESS
    @test init_cloud.stage_index == test_init_cloud.stage_index
    @test init_cloud.n_Φ         == test_init_cloud.n_Φ
    @test init_cloud.resamples   == test_init_cloud.resamples
    @test init_cloud.c           == test_init_cloud.c
    @test init_cloud.total_sampling_time == test_init_cloud.total_sampling_time
    @test init_cloud.tempering_schedule  == test_init_cloud.tempering_schedule
end
