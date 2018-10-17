using DSGE, HDF5, JLD
using BenchmarkTools
using Base.Test

path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:forecast_uncertainty_override, Nullable(false))
m <= Setting(:use_population_forecast, true)
m <= Setting(:compute_shockdec_bands, true)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")

output_vars = add_requisite_output_vars([:histpseudo, :histobs,
                                         :hist4qpseudo, :hist4qobs,
                                         :forecastpseudo, :forecastobs,
                                         :forecast4qpseudo, :forecast4qobs,
                                         :shockdecpseudo, :shockdecobs,
                                         :irfpseudo, :irfobs])

##########
# Modal
##########
fc_trial = @benchmark forecast_one(m, :mode, :none, output_vars, verbose = :none) gcsample = true
mb_trial = @benchmark compute_meansbands(m, :mode, :none, output_vars;
                                         compute_shockdec_bands = true, verbose = :none) gcsample = true
mb2m_trial = @benchmark meansbands_to_matrix(m, :mode, :none, output_vars; verbose =
                                             :none) gcsample = true

trials = [fc_trial, mb_trial, mb2m_trial]
trial_names = [:forecast_one, :compute_meansbands, :meansbands_to_matrix]
group = construct_trial_group(trials, trial_names)
group_name = "modal_meansbands"

# write_ref_trial_group(group, trial_names, group_name)
print_all_benchmarks(group, "../reference/$group_name.jld", group_name, trial_names)

#######################
# Full-distribution
#######################
@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
fc_trial = @benchmark forecast_one(m, :full, :none, output_vars, verbose = :none) gcsample = true
mb_trial = @benchmark compute_meansbands(m, :full, :none, output_vars;
                                         compute_shockdec_bands = true,
                                         verbose = :none) gcsample = true
mb2m_trial = @benchmark meansbands_to_matrix(m, :full, :none, output_vars;
                                             verbose = :none) gcsample = true

trials = [fc_trial, mb_trial, mb2m_trial]
trial_names = [:forecast_one, :compute_meansbands, :meansbands_to_matrix]
group = construct_trial_group(trials, trial_names)
group_name = "full_meansbands"

# write_ref_trial_group(group, trial_names, group_name)
print_all_benchmarks(group, "../reference/$group_name.jld", group_name, trial_names)
