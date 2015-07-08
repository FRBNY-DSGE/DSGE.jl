using Base.Test

using DSGE: M990, Gensys
include("../util.jl")


Γ
model = Model()
Γ0, Γ1, C, Ψ, Π = model.eqcond(model.Θ, model.I)
Γ1, C, impact, fmat, fwt, ywt, gev, eu, loose = gensys(complex(Γ0), complex(Γ1), C, Ψ, Π, 1 + 1e-5)



# Check output matrices against Matlab output (ε = 1e-4)
for out in ["Γ1", "C", "impact", "fmat", "fwt", "ywt", "gev"]
    if out ∈ ["Γ1", "C", "impact"]
        eval(parse("$(out)_expected = readcsv(\"models/m990/gensys/$out.csv\")"))
    else
        eval(parse("$(out)_expected = readcsv_complex(\"models/m990/gensys/$out.csv\")"))
    end
    println("### $out")
    eval(parse("@test test_matrix_eq($(out)_expected, $out)"))
end

eu_matlab = readcsv("models/m990/gensys/eu.csv")
eu_matlab = convert(Array{Int64, 1}, vec(eu_matlab))
@test eu_matlab == eu

