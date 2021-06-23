using DSGE, ModelConstructors, Dates, OrderedCollections, Test, CSV, DataFrames, Random, FileIO

regenerate_output = false

## Regime switching and full-distribution
# Initialize model objects
Random.seed!(1793 * 10)
m = Model1002("ss60"; custom_settings = Dict{Symbol, Setting}(:flexible_ait_policy_change =>
                                                              Setting(:flexible_ait_policy_change,
                                                                      false),
                                                              :add_pgap => Setting(:add_pgap, false),
                                                              :add_ygap => Setting(:add_ygap, false))) # Set to false unless you want to re-generate any saved output
df_full = CSV.read(joinpath(dirname(@__FILE__), "../reference/uncertain_altpolicy_data.csv"), DataFrame)
m <= Setting(:forecast_horizons, 20)
m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, # Have to add anticipated rates to conditional data
                                :obs_nominalrate, :obs_longrate,
                                :obs_nominalrate1, :obs_nominalrate2, :obs_nominalrate3,
                                :obs_nominalrate4, :obs_nominalrate5, :obs_nominalrate6])
m <= Setting(:cond_semi_names, [:obs_spread,
                                :obs_nominalrate, :obs_longrate,
                                :obs_nominalrate1, :obs_nominalrate2, :obs_nominalrate3,
                                :obs_nominalrate4, :obs_nominalrate5, :obs_nominalrate6])
m <= Setting(:date_forecast_start, Date(2020, 6, 30))
m <= Setting(:date_conditional_end, Date(2020, 6, 30))
m <= Setting(:forecast_ndraws, 5)
m <= Setting(:forecast_block_size, get_setting(m, :forecast_ndraws))
m <= Setting(:forecast_jstep, 1)
m <= Setting(:use_parallel_workers, false)
setup_regime_switching_inds!(m; cond_type = :full)

# For parameters, use modal ones, but also draw from prior for the COVID-19 shocks
mparas = repeat(ModelConstructors.get_values(m.parameters)', get_setting(m, :forecast_ndraws), 1)
j = length(m.parameters) # To help index for regime-switching parameters
biidc_keys = Dict()
ziid_keys = Dict()
φ_keys = Dict()
for p in m.parameters # Find the keys of the COVID-19 shocks
    if haskey(p.regimes, :value)
        for k in keys(p.regimes[:value])
            if k > 1
                global j += 1
                if p.key == :σ_φ
                    φ_keys[k] = j
                elseif p.key == :σ_biidc
                    biidc_keys[k] = j
                elseif p.key == :σ_ziid
                    ziid_keys[k] = j
                end
            end
        end
    end
end

for i in 1:size(mparas, 1)
    mparas[i, φ_keys[2]] = rand(get(m.parameters[m.keys[:σ_φ]].prior)) # draw from prior for second-regime value
    mparas[i, biidc_keys[2]] = rand(get(m.parameters[m.keys[:σ_biidc]].prior)) # draw from prior for second-regime value
    mparas[i, ziid_keys[2]] = rand(get(m.parameters[m.keys[:σ_ziid]].prior)) # draw from prior for second-regime value
    mparas[i, φ_keys[3]] = rand(get(m.parameters[m.keys[:σ_φ]].prior)) # draw from prior for third-regime value
    mparas[i, biidc_keys[3]] = rand(get(m.parameters[m.keys[:σ_biidc]].prior)) # draw from prior for third-regime value
    mparas[i, ziid_keys[3]] = rand(get(m.parameters[m.keys[:σ_ziid]].prior)) # draw from prior for third-regime value
end

# Function for adding new regimes to parameters
function covid_set_regime_vals(m::AbstractDSGEModel, n::Int)
    if n > 4
        for p in m.parameters
            if haskey(p.regimes, :value)
                for i in 5:n
                    ModelConstructors.set_regime_val!(p, i, ModelConstructors.regime_val(p, 4))
                end
            end
        end
    end
end

# Run forecast!
output_vars = [:forecastobs, :bddforecastobs]
forecast_one(m, :full, :full, output_vars, verbose = :none, params = mparas, df = df_full,
             zlb_method = :temporary_altpolicy, set_regime_vals_altpolicy = covid_set_regime_vals)

# Either test against saved output or re-generate output
output_files = get_forecast_output_files(m, :full, :full, output_vars)
if regenerate_output
    using JLD2, FileIO
    if (VERSION >= v"1.5")
        JLD2.jldopen(joinpath(dirname(@__FILE__), "../reference/automatic_tempalt_zlb_fulldist_v1p5.jld2"),
                     true, true, true, IOStream) do file
            for (k, v) in output_files
                write(file, string(k), load(v, "arr"))
            end
        end
    else
        JLD2.jldopen(joinpath(dirname(@__FILE__), "../reference/automatic_tempalt_zlb_fulldist.jld2"),
                     true, true, true, IOStream) do file
            for (k, v) in output_files
                write(file, string(k), load(v, "arr"))
            end
        end
    end
else
    refdata = (VERSION >= v"1.5") ?
        load(joinpath(dirname(@__FILE__), "../reference/automatic_tempalt_zlb_fulldist_v1p5.jld2")) :
        load(joinpath(dirname(@__FILE__), "../reference/automatic_tempalt_zlb_fulldist.jld2"))
    @testset "Automatic enforcement of ZLB as a temporary alternative policy during full-distribution forecast" begin
        for (k, v) in output_files # Test variables which don't have forward-looking measurement equations

            # measurement eqns for forward looking variables e.g. :obs_longrate are still being updated.
            if k == :bddforecastobs
                @test_matrix_approx_eq refdata[string(k)][:, # this test seems to work but perhaps due to changes in the way the automatic ZLB is
                                                                vcat(1:9, 12:13), :] load(v, "arr")[:, vcat(1:9, 12:13), :] # inferred, it no longer matches reference output. Leaving this test as "broken" for now.
                @test all(refdata[string(k)][:, m.observables[:obs_nominalrate], :] .> -1e-14)
            else
                @test @test_matrix_approx_eq refdata[string(k)][:,
                                                                vcat(1:9, 12:13), :] load(v, "arr")[:, vcat(1:9, 12:13), :]
                @test !all(refdata[string(k)][:, m.observables[:obs_nominalrate], :] .> -1e-14)
            end
        end
    end
    for v in values(output_files)
        rm(v)
    end
end
