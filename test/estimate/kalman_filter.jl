using DSGE, HDF5, Base.Test
include("../util.jl")

path = dirname(@__FILE__)

# Initialize arguments to function
h5 = h5open("$path/../reference/kalman_filter_args.h5")
for arg in ["data", "lead", "a", "F", "b", "H", "var", "z0", "vz0"]
    eval(parse("$arg = read(h5, \"$arg\")"))
end
close(h5)

lead = round(Int, lead)

m = Model990()

# Method with all arguments provided (9)
out_3 = kalman_filter(m, data, F, a, H, b, var, z0, vz0; lead = lead)
out_9 = kalman_filter(m, data, F, a, H, b, var, z0, vz0; lead = lead, allout = true)

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
out_3 = kalman_filter(m, data, F, a, H, b, var; lead = lead)
out_9 = kalman_filter(m, data, F, a, H, b, var; lead = lead, allout = true)

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
        # These matrix entries are especially large, averaging 1e5, so we allow
        # greater ϵ
        @test_matrix_approx_eq_eps expect actual 1e-1 1e-2
    else
        @test_matrix_approx_eq expect actual
    end
end
close(h5)

nothing
