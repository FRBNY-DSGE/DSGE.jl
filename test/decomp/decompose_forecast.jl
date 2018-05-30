using DSGE, Base.Test, JLD

path = dirname(@__FILE__)

# Initialize models
function make_test_model(year::Int)
    m = Model990()
    vint = Dates.format(Date(year, 4, 10), "yymmdd")
    m <= Setting(:data_vintage, vint)
    m <= Setting(:cond_vintage, vint)
    m <= Setting(:date_forecast_start, DSGE.quartertodate(year, 1))
    m <= Setting(:date_conditional_end, DSGE.quartertodate(year, 1))
    m <= Setting(:forecast_horizons, 16)
    return m
end

m_new = make_test_model(2016)
m_old = make_test_model(2014)

# Read in data and parameters
file = jldopen("$path/../reference/decompose_forecast_args.jld", "r")
df_new = read(file, "df_new")
df_old = read(file, "df_old")
params_new = read(file, "params_new")
params_old = read(file, "params_old")
close(file)

# Read in expected outputs
exp_decomps = jldopen("$path/../reference/decompose_forecast_out.jld", "r") do file
    read(file, "all_decomps")
end

cond_types = [(:none, :none), (:none, :full), (:full, :none), (:full, :full)]
all_hs = [4, 12, 1:16] # h < k, h = k, h > k, h = 1:H

for (cond_new, cond_old) in cond_types
    # Test inner method
    for hs in all_hs
        @time decomps = DSGE.decompose_forecast(m_new, m_old,
                               cond_new == :full ? df_new : df_new[1:end-1, :],
                               cond_old == :full ? df_old : df_old[1:end-1, :],
                               params_new, params_old, cond_new, cond_old,
                               [:obs, :pseudo]; hs = hs, individual_shocks = true, check = true)

        for k in keys(decomps)
            @test decomps[k] ≈ all_decomps[(cond_new, cond_old, hs)][k]
        end
    end
end

nothing