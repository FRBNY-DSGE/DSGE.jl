using DSGE, JLD, Base.Test

path = dirname(@__FILE__)

# Set up
m = AnSchorfheide(testing = true)

kal1, kal2 = jldopen("$path/../reference/kalman_cat_args.jld", "r") do file
    read(file, "kal1"), read(file, "kal2")
end

# Concatenate Kalmans
kal12 = cat(m, kal1, kal2)

# Test equality
exp_kal12 = jldopen("$path/../reference/kalman_cat_out.jld", "r") do file
    read(file, "kal12")
end

for arg in fieldnames(kal1)
    if arg == :L
        @test exp_kal12[arg] ≈ kal12[arg]
    else
        @test exp_kal12[arg] ≈ kal12[arg]
    end
end

nothing
