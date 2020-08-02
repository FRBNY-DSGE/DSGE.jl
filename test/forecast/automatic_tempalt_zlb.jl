using DSGE, ModelConstructors, Dates, OrderedCollections, Test, CSV, DataFrames, Random, FileIO

regenerate_output = false

# Test automatic enforcement of ZLB as a temporary ZLB with no regime-switching
custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_presample_start     => Setting(:date_presample_start, Date(1959, 9, 30)),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_conditional_end, DSGE.quartertodate("2016-Q3")),
    :n_mon_anticipated_shocks => Setting(:n_mon_anticipated_shocks, 6))

dfs = Dict()
path = dirname(@__FILE__)
dfs[:full] = load("$path/../reference/regime_switch_data.jld2", "full")
dfs[:none] = load("$path/../reference/regime_switch_data.jld2", "none")
dfs[:full][end, :obs_nominalrate] = .5
dfs[:none][end, :obs_nominalrate] = 0.5
dfs[:full][end, :obs_gdp] = -.5

m = Model1002("ss10", custom_settings = custom_settings)
m <= Setting(:rate_expectations_source, :ois)
m.settings[:regime_switching] = Setting(:regime_switching, false)
m.settings[:forecast_ndraws] = Setting(:forecast_ndraws, 3)
m.settings[:forecast_block_size] = Setting(:forecast_block_size, 2)
m.settings[:forecast_jstep] = Setting(:forecast_jstep, 1)
output_vars = [:forecastobs, :bddforecastobs]

para = repeat(Matrix(map(x -> x.value, m.parameters)'), 3, 1)
out = Dict()
Random.seed!(1793 * 5)
for cond_type in [:full, :none]
    out[cond_type] = Dict()
    forecast_one(m, :mode, cond_type, output_vars, verbose = :none, params = para[1, :], df = dfs[cond_type],
                 zlb_method = :temporary_altpolicy)
    forecast_one(m, :full, cond_type, output_vars, verbose = :none, params = para, df = dfs[cond_type],
                 zlb_method = :temporary_altpolicy)
    for input_type in [:mode, :full]
        out[cond_type][input_type] = Dict()
        outfn = get_forecast_output_files(m, input_type, cond_type, output_vars)
        out[cond_type][input_type][:bddforecastobs] = load(outfn[:bddforecastobs], "arr")
        out[cond_type][input_type][:forecastobs]    = load(outfn[:forecastobs], "arr")
    end
end

@testset "Automatic enforcement of ZLB as temporary alternative policy (no regime-switching)" begin
    for cond_type in [:none, :full]
        for input_type in [:mode, :full]
            if input_type == :mode
                @test all(out[cond_type][input_type][:bddforecastobs][m.observables[:obs_nominalrate], :] .> -1e-14)
            else
                @test all(out[cond_type][input_type][:bddforecastobs][:, m.observables[:obs_nominalrate], :] .> -1e-14)
            end
        end
    end
end

## Regime switching and full-distriution
# Initialize model objects
Random.seed!(1793)
m = Model1002("ss60")
if VERSION >= v"1.3"
    df_full = DataFrame!(CSV.File(joinpath(dirname(@__FILE__), "../reference/uncertain_altpolicy_data.csv")))
else
    df_full = DataFrame(CSV.read(joinpath(dirname(@__FILE__), "../reference/uncertain_altpolicy_data.csv")))
end
m <= Setting(:forecast_horizons, 30)
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
m <= Setting(:forecast_ndraws, 10)
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
    JLD2.jldopen(joinpath(dirname(@__FILE__), "../reference/automatic_tempalt_zlb_fulldist.jld2"), true, true, true, IOStream) do file
        for (k, v) in output_files
            write(file, string(k), load(v, "arr"))
        end
    end
else
    refdata = load(joinpath(dirname(@__FILE__), "../reference/automatic_tempalt_zlb_fulldist.jld2"))
    @testset "Automatic enforcement of ZLB as a temporary alternative policy during full-distribution forecast" begin
        for (k, v) in output_files
            @test @test_matrix_approx_eq refdata[string(k)] load(v, "arr")
            if k == :bddforecastobs
                @test all(refdata[string(k)][:, m.observables[:obs_nominalrate], :] .> -1e-14)
            else
                @test !all(refdata[string(k)][:, m.observables[:obs_nominalrate], :] .> -1e-14)
            end
        end
    end
    for v in values(output_files)
        rm(v)
    end
end
