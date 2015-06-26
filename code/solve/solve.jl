# Description

function solve(model::Model)
    include("gensys.jl")
    
    G0, G1, C, Ψ, Π = model.eqcond(model.Θ, model.I)

    # Currently throws a SingularException
    G1, C, impact, fmat, fwt, ywt, gev, eu, loose = GenSys.gensys(G0, G1, C, Ψ, Π, 1 + 1e-6)
end
