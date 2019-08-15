using ClusterManagers, DSGE

fn = dirname(@__FILE__)

# What do you want to do?
run_full_forecast = true
do_histforecast   = true
do_shockdecs      = true
do_irfs           = true
make_plots        = true
make_packet       = true

# Initialize model objects and desired settings
m = Model1002("ss10")
m <= Setting(:sampling_method, :SMC)
usual_settings!(m, "190808"; cdvt = "190808", fcast_date = quartertodate("2019-Q3"))
m <= Setting(:saveroot, "$(fn)/../../save/") # set save root
m <= Setting(:dataroot, "$(fn)/../../save/input_data/") # set data root
m <= Setting(:date_forecast_end, quartertodate("2030-Q1"))
forecast_string = "" # to identify different forecasts from each other
model_string = "m1002/ss10" # used further down to construct file paths

# Override input file locations for saved estimation output
overrides = forecast_input_file_overrides(m)
overrides[:full] = "$(fn)/../../save/output_data/m1002/ss10/estimate/raw/smcsave_vint=190712.h5"

# Load data from SMC cloud estimation and guarantee saved parameters are
# a draw from the weighted particle distribution
weighted_data = load("$(fn)/../../save/output_data/m1002/ss10/estimate/raw/smc_cloud_vint=190712.jld2")
Random.seed!(1793)
params = Matrix{Float64}(DSGE.get_vals(weighted_data["cloud"])[:,resample(weighted_data["W"][:,end])]')
h5open("$(fn)/../../save/output_data/m1002/ss10/estimate/raw/smcsave_vint=190712.h5", "w") do file
    write(file, "smcparams", params)
end

# Full-distribution forecast
if run_full_forecast
    if get_setting(m, :sampling_method) == :SMC
        m <= Setting(:forecast_block_size, 500)
    end

    output_vars = Vector{Symbol}(undef,0)
    if do_histforecast
        # Write data to create historical and forecast output
        output_vars = vcat(output_vars, [:histpseudo, :histobs, :histstdshocks,
                                         :hist4qpseudo, :hist4qobs, :histutpseudo,
                                         :forecastpseudo, :forecastobs, :forecastutpseudo,
                                         :forecast4qpseudo, :forecast4qobs, :forecaststdshocks])
    end
    if do_shockdecs
        # Shock decompositions of forecasts
        output_vars = vcat(output_vars, [:dettrendobs, :dettrendpseudo, :trendobs,
                                         :trendpseudo, :shockdecpseudo, :shockdecobs])
    end
    if do_irfs
        # Impulse response to all shocks, including for endogenous states
        output_vars = vcat(output_vars, [:irfstates, :irfobs, :irfpseudo])
    end

    my_procs = addprocs(20)
    @everywhere using DSGE, OrderedCollections
    usual_forecast(m, :full, :none, output_vars,
                       est_override = overrides[:full],
                       forecast_string = forecast_string,
                       density_bands = [.5, .6, .68, .7, .8, .9])
end

# Packet
if make_packet
    gr()

    output_vars = [:forecastobs, :forecastpseudo]
    if do_shockdecs
        output_vars = vcat(output_vars, [:shockdecobs, :shockdecpseudo])
    end
    if do_irfs
        # we do not plot irfs for endogenous states b/c too many!
        output_vars = vcat(output_vars, [:irfobs, :irfpseudo])
        sections = [:estimation, :forecast, :irf]
    else
        sections = [:estimation, :forecast]
    end
    if make_plots
        plot_standard_packet(m, :full, :none, output_vars,
                             forecast_string = forecast_string,
                             sections = sections)
    end

    # Choose forecast centered or the standard packet.
    # Note that these functions write to the same tex file,
    # so you can use only one for a given data vintage!
    # write_forecast_centric_packet(m, :full, :none, output_vars,
    #                               sections = sections, forecast_string = forecast_string)
    write_standard_packet(m, :full, :none, output_vars,
                                  sections = sections, forecast_string = forecast_string)
    moment_tables(m)
end

rmprocs(my_procs)
nothing
