# Description

using .AbstractModel
include("gensys.jl")

function solve(model::Model)
    G0, G1, C, Ψ, Π = model.eqcond(model.Θ, model.I)

    G1, C, impact, fmat, fwt, ywt, gev, eu, loose = GenSys.gensys(complex(G0), complex(G1), C, Ψ, Π, 1 + 1e-5)
end
