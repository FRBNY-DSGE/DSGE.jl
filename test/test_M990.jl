using Base.Test

using DSGE: DistributionsExt, AbstractModel, M990



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
    # ModelInds objecto creation
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

    # TODO: check output matrices against Matlab output
    
    println("eqcond.jl tests passed")
end
