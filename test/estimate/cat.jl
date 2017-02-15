using DSGE, JLD, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Set up
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6))
m = Model990(custom_settings = custom_settings, testing = true)

kal1, kal2 = jldopen("$path/../reference/kalman_cat_args.jld", "r") do file
    read(file, "kal1"), read(file, "kal2")
end

# Concatenate Kalmans
kal12 = cat(m, kal1, kal2; allout = true)

# Test equality
exp_kal12 = jldopen("$path/../reference/kalman_cat_out.jld", "r") do file
    read(file, "kal12")
end

for arg in fieldnames(kal1)
    if arg == :L
        @test exp_kal12[arg] ≈ kal12[arg]
    else
        @test_approx_eq exp_kal12[arg] kal12[arg]
    end
end

nothing