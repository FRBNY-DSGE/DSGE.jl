using DSGE
using HDF5, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Initialize arguments to function
h5 = h5open("$path/../reference/kalman_filter_args.h5")
for arg in ["data", "lead", "a", "F", "b", "H", "var", "z0", "vz0"]
    eval(parse("$arg = read(h5, \"$arg\")"))
end
close(h5)

lead = round(Int, lead)

# Method with all arguments provided (9)
out_3 = kalman_filter(data, lead, a, F, b, H, var, z0, vz0)
out_9 = kalman_filter(data, lead, a, F, b, H, var, z0, vz0; allout=true)

h5 = h5open("$path/../reference/kalman_filter_out9.h5")
for out in [:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse,
            :rmsd, :filt, :vfilt]
    expect = read(h5, "$out")
    actual = out_9[out]

    if out == :zend
        expect = reshape(expect, length(expect), 1)
    end

    if ndims(expect) == 0
        @test_approx_eq expect actual
    else
        @test_matrix_approx_eq expect actual
    end
end
close(h5)

# Method with optional arguments omitted (7)
out_3 = kalman_filter(data, lead, a, F, b, H, var)
out_9 = kalman_filter(data, lead, a, F, b, H, var; allout=true)

h5 = h5open("$path/../reference/kalman_filter_out7.h5")
for out in [:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse,
            :rmsd, :filt, :vfilt]
    expect = read(h5, "$out")
    actual = out_9[out]

    if out == :zend
        expect = reshape(expect, length(expect), 1)
    end

    if out == :L
        @test_approx_eq_eps expect actual 1e-4
    elseif out ∈ [:Pend, :vpred, :vfilt]
        # These matrix entries are especially large, averaging 1e5, so we allow greater ϵ
        @test_matrix_approx_eq_eps expect actual 1e-1 1e-2
    else
        @test_matrix_approx_eq expect actual
    end
end
close(h5)



# Two-part Kalman filter
data = h5open("$path/../reference/smoother_args.h5", "r") do h5
    read(h5, "data")
end

m = Model990()
m.testing = true
m <= Setting(:n_anticipated_shocks, 6)
DSGE.init_model_indices!(m)
m <= Setting(:date_forecast_start, quartertodate("2016-Q1"))

TTT, RRR, CCC = solve(m)
meas = measurement(m, TTT, RRR, CCC)
QQ, ZZ, DD = meas.QQ, meas.ZZ, meas.DD

A0 = zeros(size(TTT, 1))
P0 = QuantEcon.solve_discrete_lyapunov(TTT, RRR*QQ*RRR')

k, R1, R2, R3 = kalman_filter_2part(m, data', TTT, RRR, CCC, A0, P0; allout = true, include_presample = true)
h5 = h5open("$path/../reference/kalman_filter_2part_out.h5", "r")
for out in [:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse,
            :rmsd, :filt, :vfilt, :z0, :vz0]
    expect = read(h5, "$out")
    actual = k[out]

    @test_approx_eq expect actual
end
close(h5)



# Concatenation of Kalman objects
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

for arg in fieldnames(k)
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
    
k12 = cat(m, k1, k2; allout = true, regime_switch = false)
k23 = cat(m, k2, k3; allout = true, regime_switch = true)

for arg in fieldnames(k)
    if arg == :L
        @test d12[arg] ≈ k12[arg]
        @test d23[arg] ≈ k23[arg]
    else
        @test_approx_eq d12[arg] k12[arg]
        @test_approx_eq d23[arg] k23[arg]
    end
end


nothing
