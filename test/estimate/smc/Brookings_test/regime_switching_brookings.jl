using DSGE, OrderedCollections, CSV, DataFrames
using ClusterManagers, HDF5
# import DSGE: usual_settings!, usual_forecast!
using ModelConstructors, Dates, JLD, JLD2, SMC, StateSpaceRoutines, Nullables

include("../../../src/estimate/param_regimes.jl")
include("util_brookings.jl")

switch_all = false
switch_policy = false
switch_kappa = true
run_stored = true # Only set to true when on the Brookings_PC branch of DSGE.
## When using regime-switching parameters, set to false.

# Kalman Filter when run_stored is true
function kalman_stored(m,df,system)

data = df_to_matrix(m, df; cond_type = :none, in_sample = true)
start_date = max(date_presample_start(m), df[1, :date])
s_0 = Vector{Float64}(undef, 0)
P_0 = Matrix{Float64}(undef, 0, 0)
include_presample = true
outputs = [:loglh, :pred, :filt]
tol = 0.0

    T = size(data, 2)

    T = size(data, 2)
    if !isempty(data)
        # Calculate the number of periods since start date for each regime
        if haskey(m.settings, :n_hist_regimes)
            n_regime_periods = Vector{Int}(undef, get_setting(m, :n_hist_regimes) + get_setting(m, :n_cond_regimes))
            for k in 1:(get_setting(m, :n_hist_regimes)  + get_setting(m, :n_cond_regimes))
                n_regime_periods[k] = subtract_quarters(get_setting(m, :regime_dates)[k], start_date)
            end
        else
            n_regime_periods = Vector{Int}(undef, length(get_setting(m, :regime_dates)))
            for (k, v) in get_setting(m, :regime_dates)
                n_regime_periods[k] = subtract_quarters(v, start_date)
            end
        end

        if start_date < date_presample_start(m)
            error("Start date $start_date must be >= date_presample_start(m)")
        elseif 0 < subtract_quarters(date_zlb_start(m), start_date) < T
            regime_inds = Vector{UnitRange{Int64}}(undef, length(n_regime_periods) + 1)

            n_nozlb_periods = subtract_quarters(date_zlb_start(m), start_date) # number of periods since start date for start of ZLB

            # Get index of next regime after ZLB starts.
            # Note that it cannot be 1 b/c the first regime starts at the start date
            i_splice_zlb = findfirst(n_nozlb_periods .< n_regime_periods)

            if isnothing(i_splice_zlb)
                check_zlb_coincide, i_zlb_start = if n_nozlb_periods == n_regime_periods[end]
                    n_regime_periods[end], length(n_regime_periods)
                else
                    -1, 0
                end
            else
                i_zlb_start = i_splice_zlb - 1
                check_zlb_coincide = n_regime_periods[i_zlb_start]
            end
            if check_zlb_coincide == n_nozlb_periods
                # In this case, the post-ZLB coincides with a specific regime,
                # so we return the regime indices w/out splicing the ZLB in.
                return regime_indices(m, start_date, iterate_quarters(start_date, T - 1)), i_zlb_start, false
            else
                # Populate vector of regime indices
                if isnothing(i_splice_zlb) # post-ZLB is the last regime
                    for reg in 1:(length(n_regime_periods) - 1)
                        regime_inds[reg] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
                    end
                    regime_inds[end - 1] = (n_regime_periods[end] + 1):n_nozlb_periods
                    regime_inds[end]     = (n_nozlb_periods + 1):T

                    i_zlb_start = length(regime_inds)
                elseif i_splice_zlb > 2 # at least one full regime before ZLB starts
                    regime_inds[1] = 1:n_regime_periods[2]
                    for reg in 2:i_splice_zlb - 2 # if i_splice_zlb == 3, then this loop does not run
                        regime_inds[reg] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
                    end
                    regime_inds[i_splice_zlb - 1] = (n_regime_periods[i_splice_zlb - 1] + 1):n_nozlb_periods
                    regime_inds[i_splice_zlb]     = (n_nozlb_periods + 1):(n_regime_periods[i_splice_zlb])

                    i_zlb_start = i_splice_zlb

                    # Index reg + 1 b/c we have spliced pre- and post- ZLB regime in
                    for reg in i_splice_zlb:(length(n_regime_periods) - 1)
                        regime_inds[reg + 1] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
                    end
                    regime_inds[end] = (n_regime_periods[end] + 1):T
                else # post- ZLB regimes start in the very first regime specified by user
                    regime_inds[1] = 1:n_nozlb_periods
                    regime_inds[2] = (n_nozlb_periods + 1):n_regime_periods[2]

                    i_zlb_start = 2

                    # Index reg + 1 b/c spliced pre- and post-ZLB regime in
                    for reg in 2:(length(n_regime_periods) - 1)
                        regime_inds[reg + 1] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
                    end
                    regime_inds[end] = (n_regime_periods[end] + 1):T
                end
                splice_zlb_regime = true
            end
        else
            # This is the case that date_zlb_start <= start_date so the first regime is the post-ZLB regime (no pre-ZLB)
            # Since there is no ZLB regime to splice in, we return the regime_inds corresponding to the user-specified regimes
            regime_inds = Vector{UnitRange{Int64}}(undef, length(n_regime_periods))
            for reg in 1:(length(n_regime_periods) - 1)
                regime_inds[reg] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
            end
            regime_inds[end] = (n_regime_periods[end] + 1):T

            i_zlb_start = 1 # start immediately in post-ZLB

            splice_zlb_regime = date_zlb_start(m) == start_date ? false : true
        end
    else # Empty, so we ignore regime switching
        regime_inds = UnitRange{Int64}[1:T]
        i_zlb_start = 0
        splice_zlb_regime = false
    end
    # Remove regimes that only are for the forecast period (i.e. bigger than number of rows in dataframe)
    filter!(x -> last(x) <= T, regime_inds)

    # Partition sample into regimes, including the pre- and post-ZLB regimes.
    # regime_inds, i_zlb_start, splice_zlb_regime = zlb_plus_regime_indices(m, data, start_date)

    # Get system matrices for each regime
    n_reg = length(regime_inds)
    ind_zlb_start = i_zlb_start

    S = Float64

    if n_reg <= 0 || ind_zlb_start <= 0
        throw(DomainError())
    end

    if start_date < date_presample_start(m)
        error("Start date $start_date must be >= date_presample_start(m)")

        # TODO: This technically doesn't handle the case where the end_date of the sample
        # is before the start of the ZLB
    elseif date_presample_start(m) <= start_date <= date_zlb_start(m)
        TTTs = Vector{Matrix{S}}(undef, n_reg)
        RRRs = Vector{Matrix{S}}(undef, n_reg)
        CCCs = Vector{Vector{S}}(undef, n_reg)
        ZZs  = Vector{Matrix{S}}(undef, n_reg)
        DDs  = Vector{Vector{S}}(undef, n_reg)
        QQs  = Vector{Matrix{S}}(undef, n_reg)
        EEs  = Vector{Matrix{S}}(undef, n_reg)

        shock_inds = inds_shocks_no_ant(m)

        if splice_zlb_regime
            if ind_zlb_start == 2
                if false#n_mon_anticipated_shocks(m) == 0
                    QQ_preZLB = copy(system[1][:QQ]) # anticipated shocks are zero already
                else
                    QQ_preZLB                         = zeros(size(system[1][:QQ]))
                    QQ_preZLB[shock_inds, shock_inds] = system[1][:QQ][shock_inds, shock_inds]
                end
                QQ_postZLB = copy(system[1][:QQ])

                TTTs[1:2] = [system[1][:TTT], system[1][:TTT]]::Vector{Matrix{S}}
                RRRs[1:2] = [system[1][:RRR], system[1][:RRR]]::Vector{Matrix{S}}
                CCCs[1:2] = [system[1][:CCC], system[1][:CCC]]::Vector{Vector{S}}
                ZZs[1:2]  = [system[1][:ZZ],  system[1][:ZZ]]::Vector{Matrix{S}}
                DDs[1:2]  = [system[1][:DD],  system[1][:DD]]::Vector{Vector{S}}
                EEs[1:2]  = [system[1][:EE],  system[1][:EE]]::Vector{Matrix{S}}
                QQs[1]    = QQ_preZLB
                QQs[2]    = QQ_postZLB

                if n_reg >= 3
                    input = 2:(n_reg - 1) # minus 1 b/c spliced in a ZLB regime

                    TTTs[3:end] = [system[i][:TTT] for i in input]::Vector{Matrix{S}}
                    RRRs[3:end] = [system[i][:RRR] for i in input]::Vector{Matrix{S}}
                    CCCs[3:end] = [system[i][:CCC] for i in input]::Vector{Vector{S}}
                    ZZs[3:end]  = [system[i][:ZZ] for i in input]::Vector{Matrix{S}}
                    DDs[3:end]  = [system[i][:DD] for i in input]::Vector{Vector{S}}
                    QQs[3:end]  = [system[i][:QQ] for i in input]::Vector{Matrix{S}}
                    EEs[3:end]  = [system[i][:EE] for i in input]::Vector{Matrix{S}}
                end
            elseif ind_zlb_start >= 3
                # Populate matrices for regimes before the ZLB split
                pre_zlb_input = 1:(ind_zlb_start - 2) # minus 2 b/c ind_zlb_start - 1 will be the pre-ZLB regime
                splice_reg    = ind_zlb_start - 1

                TTTs[pre_zlb_input] = [system[i][:TTT] for i in pre_zlb_input]::Vector{Matrix{S}}
                RRRs[pre_zlb_input] = [system[i][:RRR] for i in pre_zlb_input]::Vector{Matrix{S}}
                CCCs[pre_zlb_input] = [system[i][:CCC] for i in pre_zlb_input]::Vector{Vector{S}}
                ZZs[pre_zlb_input]  = [system[i][:ZZ] for i in pre_zlb_input]::Vector{Matrix{S}}
                DDs[pre_zlb_input]  = [system[i][:DD] for i in pre_zlb_input]::Vector{Vector{S}}
                EEs[pre_zlb_input]  = [system[i][:EE] for i in pre_zlb_input]::Vector{Matrix{S}}

                if false#n_mon_anticipated_shocks(m) == 0
                    QQs[pre_zlb_input] = [system[i][:QQ] for i in pre_zlb_input]::Vector{Matrix{S}}
                else
                    for i in pre_zlb_input
                        QQ_preZLB                         = zeros(size(system[i][:QQ]))
                        QQ_preZLB[shock_inds, shock_inds] = system[i][:QQ][shock_inds, shock_inds]
                        QQs[i]                            = QQ_preZLB
                    end
                end

                # Regime in which the ZLB occurs
                if false#n_mon_anticipated_shocks(m) == 0
                    QQ_preZLB = copy(system[splice_reg, :QQ]) # anticipated shocks are zero already
                else
                    QQ_preZLB                         = zeros(size(system[splice_reg][:QQ]))
                    QQ_preZLB[shock_inds, shock_inds] = system[splice_reg][:QQ][shock_inds, shock_inds]
                end
                QQ_postZLB = copy(system[splice_reg][:QQ])
                TTTs[splice_reg:ind_zlb_start] = [system[splice_reg][:TTT], system[splice_reg][:TTT]]::Vector{Matrix{S}}
                RRRs[splice_reg:ind_zlb_start] = [system[splice_reg][:RRR], system[splice_reg][:RRR]]::Vector{Matrix{S}}
                CCCs[splice_reg:ind_zlb_start] = [system[splice_reg][:CCC], system[splice_reg][:CCC]]::Vector{Vector{S}}
                ZZs[splice_reg:ind_zlb_start]  = [system[splice_reg][:ZZ],  system[splice_reg][:ZZ]]::Vector{Matrix{S}}
                DDs[splice_reg:ind_zlb_start]  = [system[splice_reg][:DD],  system[splice_reg][:DD]]::Vector{Vector{S}}
                EEs[splice_reg:ind_zlb_start]  = [system[splice_reg][:EE],  system[splice_reg][:EE]]::Vector{Matrix{S}}
                QQs[splice_reg]     = QQ_preZLB
                QQs[ind_zlb_start] = QQ_postZLB

                # Handle regimes after the ZLB split
                if n_reg >= ind_zlb_start + 1
                    post_zlb_input = ind_zlb_start:(n_reg - 1) # minus 1 to both b/c spliced in the ZLB regime

                    TTTs[4:end] = [system[i][:TTT] for i in post_zlb_input]::Vector{Matrix{S}}
                    RRRs[4:end] = [system[i][:RRR] for i in post_zlb_input]::Vector{Matrix{S}}
                    CCCs[4:end] = [system[i][:CCC] for i in post_zlb_input]::Vector{Vector{S}}
                    ZZs[4:end]  = [system[i][:ZZ] for i in post_zlb_input]::Vector{Matrix{S}}
                    DDs[4:end]  = [system[i][:DD] for i in post_zlb_input]::Vector{Vector{S}}
                    QQs[4:end]  = [system[i][:QQ] for i in post_zlb_input]::Vector{Matrix{S}}
                    EEs[4:end]  = [system[i][:EE] for i in post_zlb_input]::Vector{Matrix{S}}
                end
             end
        else
            if ind_zlb_start > 1
                pre_zlb_input = 1:(ind_zlb_start - 1) # minus 1 b/c post-ZLB coincides with a regime

                TTTs[pre_zlb_input] = [system[i][:TTT] for i in pre_zlb_input]::Vector{Matrix{S}}
                RRRs[pre_zlb_input] = [system[i][:RRR] for i in pre_zlb_input]::Vector{Matrix{S}}
                CCCs[pre_zlb_input] = [system[i][:CCC] for i in pre_zlb_input]::Vector{Vector{S}}
                ZZs[pre_zlb_input]  = [system[i][:ZZ] for i in pre_zlb_input]::Vector{Matrix{S}}
                DDs[pre_zlb_input]  = [system[i][:DD] for i in pre_zlb_input]::Vector{Vector{S}}
                EEs[pre_zlb_input]  = [system[i][:EE] for i in pre_zlb_input]::Vector{Matrix{S}}

                if false#n_mon_anticipated_shocks(m) == 0
                    QQs[pre_zlb_input] = Matrix{S}[system[i][:QQ] for i in pre_zlb_input]::Vector{Matrix{S}}
                else
                    for i in pre_zlb_input
                        QQ_preZLB                         = zeros(size(system[i][:QQ]))
                        QQ_preZLB[shock_inds, shock_inds] = system[i][:QQ][shock_inds, shock_inds]
                        QQs[i]                            = QQ_preZLB
                    end
                end
            end

            # Handle regimes after the ZLB split
            post_zlb_input = ind_zlb_start:n_reg # no minus 1 b/c didn't splice in the regime

            TTTs[post_zlb_input] = [system[i][:TTT] for i in post_zlb_input]::Vector{Matrix{S}}
            RRRs[post_zlb_input] = [system[i][:RRR] for i in post_zlb_input]::Vector{Matrix{S}}
            CCCs[post_zlb_input] = [system[i][:CCC] for i in post_zlb_input]::Vector{Vector{S}}
            ZZs[post_zlb_input]  = [system[i][:ZZ] for i in post_zlb_input]::Vector{Matrix{S}}
            DDs[post_zlb_input]  = [system[i][:DD] for i in post_zlb_input]::Vector{Vector{S}}
            QQs[post_zlb_input]  = [system[i][:QQ] for i in post_zlb_input]::Vector{Matrix{S}}
            EEs[post_zlb_input]  = [system[i][:EE] for i in post_zlb_input]::Vector{Matrix{S}}
         end
    elseif date_zlb_start(m) < start_date
        # We always use post-ZLB matrices in this case,
        # hence we just return the matrices of the `RegimeSwitchingSystem`
        input = 1:n_reg
        TTTs = [system[i][:TTT] for i in input]::Vector{Matrix{S}}
        RRRs = [system[i][:RRR] for i in input]::Vector{Matrix{S}}
        CCCs = [system[i][:CCC] for i in input]::Vector{Vector{S}}
        ZZs  = [system[i][:ZZ] for i in input]::Vector{Matrix{S}}
        DDs  = [system[i][:DD] for i in input]::Vector{Vector{S}}
        QQs  = [system[i][:QQ] for i in input]::Vector{Matrix{S}}
        EEs  = [system[i][:EE] for i in input]::Vector{Matrix{S}}
    end
#=
    TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs = zlb_plus_regime_matrices(m, system, length(regime_inds),
                                                                    start_date;
                                                                    ind_zlb_start = i_zlb_start,
                                                                    splice_zlb_regime = splice_zlb_regime)
=#
    # If s_0 and P_0 provided, check that rows and columns corresponding to
    # anticipated shocks are zero in P_0
    if !isempty(s_0) && !isempty(P_0)
        ant_state_inds = setdiff(1:n_states_augmented(m), inds_states_no_ant(m))
        @assert all(x -> x == 0, P_0[:, ant_state_inds])
        @assert all(x -> x == 0, P_0[ant_state_inds, :])
    end

    # Specify number of presample periods if we don't want to include them in
    # the final results
    Nt0 = include_presample ? 0 : n_presample_periods(m)

    # Run Kalman filter, construct Kalman object, and return
    out = kalman_filter(regime_inds, data, TTTs, RRRs, CCCs, QQs,
                        ZZs, DDs, EEs, s_0, P_0; outputs = outputs,
                        Nt0 = Nt0, tol = tol)
return out[1]
end


# Initialize model objects
function model_init(;subspec="ss10", run_stored = false)
m = Model1002(subspec)
fp = dirname(@__FILE__)

if subspec == "ss22" || subspec == "ss24"
    standard_spec!(m, "191118", fcast_date = Date(2019, 12, 31))
    m <= Setting(:add_laborshare_measurement, false)
    m <= Setting(:hours_first_observable, false)

    m <= Setting(:date_regime2_start, Date(1990, 3, 31))
    m <= Setting(:date_regime2_start_text, "900331", true, "reg2start", "The text version to be saved of when regime 2 starts")
end

m <= Setting(:data_vintage, "191118")
m <= Setting(:saveroot, joinpath(fp, "output_data/"))
m <= Setting(:dataroot, joinpath(fp, "input_data/"))
m <= Setting(:date_forecast_start, Date(2019,12,31))

df = load_data(m, check_empty_columns = false)

m <= Setting(:n_particles, 500)#15000)
m <= Setting(:n_smc_blocks, 5)
m <= Setting(:n_mh_steps_smc, 1)
m <= Setting(:use_parallel_workers, false)
m <= Setting(:step_size_smc, 0.5)
m <= Setting(:resampler_smc, :systematic)
m <= Setting(:target_accept, 0.25)
m <= Setting(:mixture_proportion, 0.9)
m <= Setting(:use_fixed_schedule, false)
m <= Setting(:adaptive_tempering_target_smc, 0.97)
m <= Setting(:resampling_threshold, .5)
m <= Setting(:use_fixed_schedule, false)

# df_first = df[Date(1964, 3, 1) .<= df[:date]  .<= Date(1989, 12, 31), :]
# df_second = df[Date(1990, 3, 1) .<= df[:date], :]


# my_procs = addprocs_frbny(200)
# @everywhere using DSGE, OrderedCollections, DSGEModels, ModelConstructors
# @everywhere include("../../../src/estimate/param_regimes.jl")

# COMBINED
if run_stored
    m <= Setting(:period, "combined", true, "period", "period for this exercse (first or second)")
end
m <= Setting(:date_presample_start, DSGE.quartertodate("1964-Q1")) # need to reset presample
m <= Setting(:date_mainsample_start, DSGE.quartertodate("1964-Q3"))

# Set model regimes
# m <= Setting(:regime_switching, true)
# m <= Setting(:n_regimes, 2)
if !run_stored || subspec in ["ss22", "ss24"]
    m <= Setting(:regime_switching, true)
    m <= Setting(:n_regimes, 2)
    m <= Setting(:regime_dates, Dict(1 => Date(1964,3,31), 2 => Date(1990, 3, 31)))
end
if !run_stored
    setup_regime_switching_inds!(m)
end
return m, df
end

if switch_all

if run_stored
    m, df = model_init(run_stored = run_stored)

    cloud10_reg1 = JLD.load("m1002/ss10/estimate/raw/smc_cloud_npart=15000_period=r1_preZLB=false_reg2start=900331_vint=191118.jld2", "cloud").particles
    cloud10_reg2 = JLD.load("m1002/ss10/estimate/raw/smc_cloud_npart=15000_period=r2_preZLB=false_reg2start=900331_vint=191118.jld2", "cloud").particles

    inds = argmax(SMC.get_loglh(cloud10_reg1))
    para1 = cloud10_reg1[inds, 1:SMC.ind_para_end(size(cloud10_reg1, 2))]
    para2 = cloud10_reg2[inds, 1:SMC.ind_para_end(size(cloud10_reg2, 2))]

    for para in 1:length(m.parameters)
        m.parameters[para].value = para1[para]
    end

    df = df[Date(1964, 3, 1) .<= df[:date], :]
    df_first = df[Date(1964, 3, 1) .<= df[:date]  .<= Date(1989, 12, 31), :]
    df_second = df[Date(1990, 3, 1) .<= df[:date], :]
    system = compute_system(m)
    kal = DSGE.filter(m, df_first, system)
    loglh1 = kal[:loglh]
    # loglh1 = DSGE.forecast_one_draw(m, :mode, :none, [:forecastobs], para1, df_first,
    #                       return_loglh = true)

    old_system1 = compute_system(m)

    m <= Setting(:period, "second", true, "period", "period for this exercse (first or second)")
    m <= Setting(:date_presample_start, DSGE.quartertodate("1990-Q1")) # need to reset presample
    m <= Setting(:date_mainsample_start, DSGE.quartertodate("1990-Q3"))

    for para in 1:length(m.parameters)
        m.parameters[para].value = para2[para]
    end
    eq_old2 = eqcond(m)
    old_system2_1 = compute_system(m)

    m2_1 = deepcopy(m)

    system = compute_system(m)
    kal = DSGE.filter(m, df_second, system)
    loglh2 = kal[:loglh]

    # loglh2 = DSGE.forecast_one_draw(m, :mode, :none, [:forecastobs], para2, df_second,
    #                       return_loglh = true)
    m2 = deepcopy(m)
    old_system2 = compute_system(m)
end
m, df = model_init(subspec = "ss10")
# Add map from model regimes to parameter regimes
param_mat = repeat([1 2], length(m.parameters))
setup_param_regimes!(m, param_mat = param_mat)

# Set parameter values equal to Brookings cloud
cloud10_reg1 = JLD.load("m1002/ss10/estimate/raw/smc_cloud_npart=15000_period=r1_preZLB=false_reg2start=900331_vint=191118.jld2", "cloud").particles
cloud10_reg2 = JLD.load("m1002/ss10/estimate/raw/smc_cloud_npart=15000_period=r2_preZLB=false_reg2start=900331_vint=191118.jld2", "cloud").particles

inds = argmax(SMC.get_loglh(cloud10_reg1))
para1 = cloud10_reg1[inds, 1:SMC.ind_para_end(size(cloud10_reg1, 2))]
para2 = cloud10_reg2[inds, 1:SMC.ind_para_end(size(cloud10_reg2, 2))]

if !run_stored
    loglh1 = SMC.get_loglh(cloud10_reg1)[inds]
    loglh2 = SMC.get_loglh(cloud10_reg2)[inds]
end
# @show loglh1, loglh2

# create parameter regime values (initialized as
## being the same in both regimes)
paramed = copy(para1)
for p2 in 1:length(m.parameters)
    p = m.parameters[p2]
    for i in 1:2
        ModelConstructors.set_regime_prior!(p, i, p.prior)
        ModelConstructors.set_regime_fixed!(p, i, p.fixed)
        ModelConstructors.set_regime_valuebounds!(p, i, p.valuebounds)
    end
    ModelConstructors.set_regime_val!(p, 1, para1[p2])
    ModelConstructors.set_regime_val!(p, 2, para2[p2])

    append!(paramed, para2[p2])
end

new_system_1 = compute_system(m)
m_new = deepcopy(m)
# Call function forecast_one_draw with return_loglh = true
kal = DSGE.filter(m, df, new_system_1)
loglh = kal[:loglh]
# loglh = DSGE.forecast_one_draw(m, :mode, :none, [:forecastobs], paramed, df,
#                           return_loglh = true)

new_system = compute_system(m)

# Compare this value with loglh from cloud for that parameter
@show sum(loglh[1:104]) #≈ loglh1
@show sum(loglh[105:end]) #≈ loglh2: Earlier Skipping 1990-Q1 and Q2 as before mainsample
end

if switch_policy
#################################
# Only policy rule param change #
#################################
# Set parameter values equal to Brookings cloud
cloud24_reg = JLD.load("m1002/ss24/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_incZLB_reg2start=900331_reg=2_vint=191118.jld2", "cloud").particles

inds = argmax(SMC.get_loglh(cloud24_reg))
para1 = cloud24_reg[inds, 1:SMC.ind_para_end(size(cloud24_reg, 2))]

if run_stored
m, df = model_init(subspec = "ss24", run_stored = true)

for para in 1:length(m.parameters)
    m.parameters[para].value = para1[para]
end

system = compute_system(m, regime_switching = true, n_regimes = 2)
loglh1 = kalman_stored(m,df,system)

#kal = filter(m, df, system)
#loglh1 = kal[:loglh]

JLD2.jldopen("output_data/loglh_policy.jld2", "w") do file
    write(file, "loglh", loglh1)  # alternatively, say "@write file A"
    write(file, "para", m.parameters)
end
else
    m, df = model_init(subspec = "ss10", run_stored = false)
    m <= Setting(:use_population_forecast, true)
#=
    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed=false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label="\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label="\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed=false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label="\\psi_3")

    m <= parameter(:ρ, 0.7126, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.75, 0.10), fixed=false,
                   description="ρ: The degree of inertia in the monetary policy rule.",
                   tex_label="\\rho_R")
=#
# create parameter regime values (initialized as
## being the same in both regimes)
Ψ_ind = 15#(1:length(m.parameters))[[m.parameters[i].key == :ψ1 for i in 1:length(m.parameters)]]
Ψ2_ind = 16
ψ3_ind = 17
ρ_ind = 20
ρ2_ind = 21
Ψ_r2_ind = 22
Ψ2_r2_ind = 23
Ψ3_r2_ind = 24

# Add map from model regimes to parameter regimes
param_mat = repeat([1 1], length(m.parameters))
for k in vcat(15:17,20)
    param_mat[k,:] = [1 2]
end
setup_param_regimes!(m, param_mat = param_mat)


# para2 = cloud10_reg2[inds, 1:SMC.ind_para_end(size(cloud10_reg2, 2))]

# loglh1 = SMC.get_loglh(cloud24_reg)[inds]
# loglh2 = SMC.get_loglh(cloud10_reg2)[inds]

# @show loglh1, loglh2

params = para1[vcat(1:20, 25:length(para1))]
for p2 in 1:length(para1[vcat(1:20, 25:length(para1))])
    p = m.parameters[p2]
    end_reg = p2 in [15,16,17,20] ? 2 : 1
    for i in 1:end_reg
        ModelConstructors.set_regime_prior!(p, i, p.prior)
        ModelConstructors.set_regime_fixed!(p, i, p.fixed)
        ModelConstructors.set_regime_valuebounds!(p, i, p.valuebounds)
    end
    if p2 > 20
        ModelConstructors.set_regime_val!(p, 1, para1[p2+4])
    else
        ModelConstructors.set_regime_val!(p, 1, para1[p2])
    end
    if end_reg == 1
        continue
    elseif p2 < 20
        ModelConstructors.set_regime_val!(p, 2, para1[p2+7])
        append!(params, para1[p2+7])
    else
        ModelConstructors.set_regime_val!(p, 2, para1[p2+1])
        append!(params, para1[p2+1])
    end
end

# Call function forecast_one_draw with return_loglh = true
system = compute_system(m)
kal = DSGE.filter(m, df, system)
loglh = kal[:loglh]
#=
loglh = DSGE.forecast_one_draw(m, :mode, :none, [:forecastobs], params, df,
                          return_loglh = true)
=#
# Compare this value with loglh from cloud for that parameter
@show sum(loglh)
# @show sum(loglh[1:104]) #≈ loglh1
# @show sum(loglh[105:end]) #≈ loglh2
end
end

if switch_kappa
#################################
### Only kappa_p param change ###
#################################
# Set parameter values equal to Brookings cloud
cloud2_reg = JLD.load("m1002/ss22/estimate/raw/smc_cloud_friday=true_npart=20000_period=full_incZLB_reg2start=900331_reg=2_vint=191118.jld2", "cloud").particles

inds = argmax(SMC.get_loglh(cloud2_reg))
para1 = cloud2_reg[inds, 1:SMC.ind_para_end(size(cloud2_reg, 2))]

if run_stored
m, df = model_init(subspec = "ss22", run_stored = true)

for para in 1:length(m.parameters)
    m.parameters[para].value = para1[para]
end

system = compute_system(m, regime_switching = true, n_regimes = 2)
loglh1 = kalman_stored(m, df, system)
# kal = DSGE.filter(m, df, system)
# loglh1 = kal[:loglh]

JLD2.jldopen("output_data/loglh_kappa.jld2", "w") do file
    write(file, "loglh", loglh1)  # alternatively, say "@write file A"
    write(file, "para", m.parameters)
end
else


m, df = model_init(subspec = "ss10")

m <= parameter(:ζ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ζ_p: The Calvo parameter (Regime 1). In every period, intermediate goods producers optimize prices with probability (1-ζ_p). With probability ζ_p, prices are adjusted according to a weighted average of the previous period's inflation (π_t1) and steady-state inflation (π_star).",
                   tex_label="\\zeta_p")

# create parameter regime values (initialized as
## being the same in both regimes)
ζ_p_ind = 2

# Add map from model regimes to parameter regimes
param_mat = repeat([1 1], length(m.parameters))
for k in [2]
    param_mat[k,:] = [1 2]
end
setup_param_regimes!(m, param_mat = param_mat)


# para2 = cloud10_reg2[inds, 1:SMC.ind_para_end(size(cloud10_reg2, 2))]

# loglh1 = SMC.get_loglh(cloud2_reg)[inds]
# loglh2 = SMC.get_loglh(cloud10_reg2)[inds]

# @show loglh1, loglh2

params = para1[vcat(1:20, 22:length(para1))]
for p2 in 1:length(para1[vcat(1:20, 22:length(para1))])
    p = m.parameters[p2]
    end_reg = p2 in [2] ? 2 : 1
    for i in 1:end_reg
        ModelConstructors.set_regime_prior!(p, i, p.prior)
        ModelConstructors.set_regime_fixed!(p, i, p.fixed)
        ModelConstructors.set_regime_valuebounds!(p, i, p.valuebounds)
    end
    if p2 > 20
        ModelConstructors.set_regime_val!(p, 1, para1[p2+1])
    else
        ModelConstructors.set_regime_val!(p, 1, para1[p2])
    end
    if end_reg == 1
        continue
    else
        ModelConstructors.set_regime_val!(p, 2, para1[21])
        append!(params, para1[21])
    end
end

# Call function forecast_one_draw with return_loglh = true
system = compute_system(m)
kal = DSGE.filter(m, df, system)
loglh_kappa = kal[:loglh]
#=
loglh = DSGE.forecast_one_draw(m, :mode, :none, [:forecastobs], params, df,
                          return_loglh = true)
=#
# Compare this value with loglh from cloud for that parameter
@show sum(loglh_kappa)
# @show sum(loglh[1:104]) #≈ loglh1
# @show sum(loglh[105:end]) #≈ loglh2
end
end

#=
@show rawpath(m, "estimate", "smc_cloud", ["new_brookings"])
smc2(m, df_to_matrix(m, df), filestring_addl = ["new_brookings_fixed_sched"], regime_switching = true)

# FIRST
m <= Setting(:period, "first", true, "period", "period for this exercse (first or second)")
m <= Setting(:date_presample_start, DSGE.quartertodate("1964-Q1")) # need to reset presample
m <= Setting(:date_mainsample_start, DSGE.quartertodate("1964-Q3"))
smc2(m, df_to_matrix(m, df_first), filestring_addl = ["old_brookings_first"])

# SECOND
m <= Setting(:period, "second", true, "period", "period for this exercse (first or second)")
m <= Setting(:date_presample_start, DSGE.quartertodate("1990-Q1")) # need to reset presample
m <= Setting(:date_mainsample_start, DSGE.quartertodate("1990-Q3"))
smc2(m, df_to_matrix(m, df_second), filestring_addl = ["old_brookings_second"])
=#
# rmprocs(my_procs)
