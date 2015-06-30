using Base.Test

using DSGE: DistributionsExt, AbstractModel, M990
include("../util.jl")


# src/models/M990.jl
function test_all()
    test_parameters()
    test_modelinds()
    test_eqcond()

    # Model instantiation
    model = Model990()
    @test isa(model, Model)
    
    println("All tests in M990.jl passed")
end



# src/models/m990/parameters.jl
function test_parameters()
    # Parameters990 object creation
    Θ = Parameters990()
    @test isa(Θ, Parameters990)

    # TODO: Check that parameters match para, bounds, etc. vectors from Matlab
    
    println("parameters.jl tests passed")
end



# src/models/m990/modelinds.jl
function test_modelinds()
    # ModelInds object creation
    I = ModelInds990()
    @test isa(I, ModelInds)

    # Endogenous states
    endo = I.endostates
    @test length(endo) == 66
    @test endo["E_z"] == 60

    # Exogenous shocks
    exo = I.exoshocks
    @test length(exo) == 22
    @test exo["pce_sh"] == 16

    # Expectation shocks
    ex = I.expshocks
    @test length(ex) == 13
    @test ex["Erk_f_sh"] == 13

    # Equations
    eq = I.equations
    @test length(eq) == 66
    @test eq["eq_Ez"] == 60
    
    println("modelinds.jl tests passed")
end




# src/models/m990/eqcond.jl
function test_eqcond()
    # eqcond function executes successfully
    Θ = Parameters990()
    I = ModelInds990()
    G0, G1, C, Ψ, Π = eqcond(Θ, I)

    # Matrices are of expected dimensions
    @test size(G0) == (66, 66)
    @test size(G1) == (66, 66)
    @test size(C) == (66, 1)
    @test size(Ψ) == (66, 22)
    @test size(Π) == (66, 13)

    # Check output matrices against Matlab output (ε = 1e-4)
    G0_matlab = readcsv("m990/G0.csv")
    println("### G0")
    @test test_matrix_eq(G0_matlab, G0)

    G1_matlab = readcsv("m990/G1.csv")
    println("### G1")
    @test test_matrix_eq(G1_matlab, G1)

    C_matlab = readcsv("m990/C.csv")
    println("### C")
    @test test_matrix_eq(C_matlab, C)
    
    Ψ_matlab = readcsv("m990/PSI.csv")
    println("### Ψ")
    @test test_matrix_eq(Ψ_matlab, Ψ)

    Π_matlab = readcsv("m990/PIE.csv")
    println("### Π")
    @test test_matrix_eq(Π_matlab, Π)

    println("eqcond.jl tests passed")
end
