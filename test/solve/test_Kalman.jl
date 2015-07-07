using Base.Test

using DSGE.Kalman
include("../util.jl")



# Initialize arguments
# These come from the first call to the Kalman filter in gibb (line 37) -> objfcndsge (15)
#   -> dsgelh (136) -> kalcvf2NaN
for arg in ["data", "F", "b", "H", "var", "z0", "vz0"]
    eval(parse("$arg = readcsv(\"kalcvf2NaN/args/$arg.csv\")"))
end
lead = 1
a = zeros(size(F, 2), 1)


# Method with all arguments provided (9)
L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt =
    kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0)

for out in ["L", "zend", "Pend", "pred", "vpred", "yprederror", "ystdprederror", "rmse",
            "rmsd", "filt", "vfilt" ]
    eval(parse("$(out)_expected = readcsv(\"kalcvf2NaN/out9/$out.csv\")"))

    # vpred and vfilt are 3D arrays, which aren't written to csv nicely by Matlab
    if out ∈ ["vpred", "vfilt"]
        eval(parse("$(out)_expected = reshape($(out)_expected, 72, 72, 2)"))
    end

    println("### $out")
    eval(parse("@test test_matrix_eq($(out)_expected, $out)"))
end


# Method with optional arguments omitted (7)
L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt =
    kalcvf2NaN(data, lead, a, F, b, H, var)

for out in ["L", "zend", "Pend", "pred", "vpred", "yprederror", "ystdprederror", "rmse",
            "rmsd", "filt", "vfilt" ]
    eval(parse("$(out)_expected = readcsv(\"kalcvf2NaN/out7/$out.csv\")"))

    if out ∈ ["vpred", "vfilt"]
        eval(parse("$(out)_expected = reshape($(out)_expected, 72, 72, 2)"))
    end

    println("### $out")
    if out ∈ ["Pend", "vpred", "vfilt"]
        # These matrix entries are especially large, averaging 1e5, so we allow greater ε
        eval(parse("@test test_matrix_eq($(out)_expected, $out, 0.1)"))
    else
        eval(parse("@test test_matrix_eq($(out)_expected, $out)"))
    end
end


println("### kalcvf2NaN tests passed\n")
