using DSGE, ModelConstructors, HDF5, Random, JLD2, FileIO, SMC, Test, Distributions, StateSpaceRoutines, Dates

path = dirname(@__FILE__)
writing_output = false

if VERSION < v"1.5"
    ver = "111"
else
    ver = "150"
end

regenerate_sim_data = false
#########################
# Regime-switching test #
# w/ An-Schorfheide     #
#########################
m = AnSchorfheide()

# Set up SMC
m <= Setting(:sampling_method, :SMC)
m <= Setting(:data_vintage, "210101")
m <= Setting(:cond_vintage, "210101")
m <= Setting(:n_Φ, 100)
m <= Setting(:λ, 3.0)
m <= Setting(:n_particles, 500)
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

sim_filepath = joinpath(path, "..", "..", "reference/sim_data_regswitch_anschorfheide.h5")
if regenerate_sim_data
    simulate_data(m; filepath = sim_filepath)
end
data = h5read(sim_filepath, "data")

true_lik = DSGE.likelihood(m, data)

m <= Setting(:saveroot, normpath(joinpath(dirname(@__FILE__), "save")))

# With regime-switching
m <= Setting(:regime_switching, true, true, "rs2", "") # For file output purposes
m <= Setting(:n_regimes, 3)
m <= Setting(:regime_dates, Dict(1 => date_presample_start(m), 2 => DSGE.iterate_quarters(Date(1960, 3, 31), 50),
                                 3 => DSGE.iterate_quarters(Date(1960, 3, 31), 100)))
setup_regime_switching_inds!(m)
m <= Setting(:model2para_regimes, Dict{Int, Int}(1 => 1, 2 => 1, 3 => 2))

for i in 1:length(m.parameters)
    for k in 1:2
        ModelConstructors.set_regime_val!(m.parameters[i], k, m.parameters[i].value)
    end
end

param_mat = repeat([1 1 2], length(m.parameters))

DSGE.setup_param_regimes!(m, param_mat)

sys = compute_system(m)

regswitch_lik = DSGE.likelihood(m, data)

@testset "Estimate regime-switching AnSchorfheide with SMC (approx. 3 min.)" begin

    @test true_lik == regswitch_lik

    # Regime switching estimation
    true_para = ModelConstructors.get_values(m.parameters)
    Random.seed!(1793)
    DSGE.smc2(m, data, regime_switching = true, run_csminwel = false,
              verbose = :none)

    posterior_means = vec(mean(load_draws(m, :full), dims = 1))

    @test length(posterior_means) == length(m.parameters) * 2
    inds = vcat(1:4, 7:16, 18:20, 24:27, 29:length(posterior_means))
    @test maximum(abs.(true_para[inds] - posterior_means[inds])) < .6
    inds2 = [5, 17, 21, 23]
    @test maximum(abs.(true_para[inds] - posterior_means[inds])) < 1.4
end
