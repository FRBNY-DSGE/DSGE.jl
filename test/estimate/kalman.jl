import DSGE: kalcvf2NaN
using HDF5, Base.Test
include("../util.jl")

path = dirname(@__FILE__)


# Initialize arguments
# These come from the first call to the Kalman filter in gibb (line 37) -> objfcndsge (15)
#   -> dsgelh (136) -> kalcvf2NaN
# See kalcvf2NaN/test_kalcvf2NaN.m


h5 = h5open("$path/../reference/kalcvf2NaN_args.h5")
for arg in ["data", "lead", "a", "F", "b", "H", "var", "z0", "vz0"]
    eval(parse("$arg = read(h5, \"$arg\")"))
end
for arg in ["a", "b", "z0"]
    eval(parse("$arg = reshape(read(h5, \"$arg\"), length($arg), 1)"))
end

lead = round(Int,lead)


# Method with all arguments provided (9)
out_3 = kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0)
out_9 = kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0; allout=true)


#TODO clean up to use less eval
h5 = h5open("$path/../reference/kalcvf2NaN_out9.h5")
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
out_3 = kalcvf2NaN(data, lead, a, F, b, H, var)
out_9 = kalcvf2NaN(data, lead, a, F, b, H, var; allout=true)


h5 = h5open("$path/../reference/kalcvf2NaN_out7.h5")
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

nothing
