using DSGE, HDF5, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Set up 
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6))
m = Model990(custom_settings = custom_settings)
m.testing = true

d1 = Dict{Symbol, Any}()
d2 = Dict{Symbol, Any}()
d3 = Dict{Symbol, Any}()
d12 = Dict{Symbol, Any}()
d23 = Dict{Symbol, Any}()

h1 = h5open("$path/../reference/kalman_cat_args1.h5", "r")
h2 = h5open("$path/../reference/kalman_cat_args2.h5", "r")
h3 = h5open("$path/../reference/kalman_cat_args3.h5", "r")
h12 = h5open("$path/../reference/kalman_cat_out12.h5", "r")
h23 = h5open("$path/../reference/kalman_cat_out23.h5", "r")

for arg in map(symbol, names(h1))
    d1[arg] = read(h1, "$arg")
    d2[arg] = read(h2, "$arg")
    d3[arg] = read(h3, "$arg")
    d12[arg] = read(h12, "$arg")
    d23[arg] = read(h23, "$arg")
end

close(h1)
close(h2)
close(h3)
close(h12)
close(h23)

k1 = DSGE.Kalman(d1[:L], d1[:zend], d1[:Pend], d1[:pred], d1[:vpred],
    d1[:yprederror], d1[:ystdprederror], d1[:rmse], d1[:rmsd], d1[:filt],
    d1[:vfilt], d1[:z0], d1[:vz0])
k2 = DSGE.Kalman(d2[:L], d2[:zend], d2[:Pend], d2[:pred], d2[:vpred],
    d2[:yprederror], d2[:ystdprederror], d2[:rmse], d2[:rmsd], d2[:filt],
    d2[:vfilt], d2[:z0], d2[:vz0])
k3 = DSGE.Kalman(d3[:L], d3[:zend], d3[:Pend], d3[:pred], d3[:vpred],
    d3[:yprederror], d3[:ystdprederror], d3[:rmse], d3[:rmsd], d3[:filt],
    d3[:vfilt], d3[:z0], d3[:vz0])

# Concatenate Kalmans
k12 = cat(m, k1, k2; allout = true, regime_switch = false)
k23 = cat(m, k2, k3; allout = true, regime_switch = true)

# Test equality
for arg in fieldnames(k1)
    if arg == :L
        @test d12[arg] ≈ k12[arg]
        @test d23[arg] ≈ k23[arg]
    else
        @test_approx_eq d12[arg] k12[arg]
        @test_approx_eq d23[arg] k23[arg]
    end
end

nothing