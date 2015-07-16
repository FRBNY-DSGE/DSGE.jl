using Base.Test, Distributions

using DSGE: DistributionsExt, AbstractModel, M990
include("../util.jl")


# src/models/M990.jl
function test_all()
    test_parameters()
    test_modelinds()
    test_eqcond()
    test_measurement()

    # Model instantiation
    model = Model()
    @test isa(model, Model)
    
    println("### All tests in M990.jl passed\n")
end



# src/models/m990/parameters.jl
function test_parameters()
    # Parameters990 object creation
    Θ = Parameters990(spec_vars())
    @test isa(Θ, Parameters990)

    # Parameters match para, bounds, etc. vectors from Matlab (ε = 1e-4)
    para = zeros(82, 1)
    bounds = zeros(82, 2)
    pshape = zeros(82, 1)
    pmean = zeros(82, 1)
    pstdd = zeros(82, 1)
    trspec = zeros(82, 4)
    
    ignore = [Θ.del, Θ.law, Θ.epsp, Θ.epsw, Θ.gstar] # not all Params appear in para vector
    i = 1
    for φ = Θ
        if φ ∈ ignore
            continue
        end

        para[i] = φ.value
        
        (left, right) = φ.bounds
        bounds[i, 1] = left
        bounds[i, 2] = right

        if isa(φ.prior, Distributions.InverseGamma)
            pshape[i] = 4
            (α, β) = params(φ.prior)
            ν = 2α
            σ = sqrt(β/α)
            pmean[i] = σ
            pstdd[i] = ν
        else
            if isa(φ.prior, Distributions.Beta)
                pshape[i] = 1
            elseif isa(φ.prior, Distributions.Gamma)
                pshape[i] = 2
            elseif isa(φ.prior, Distributions.Normal)
                pshape[i] = 3
            end
            pmean[i] = mean(φ.prior)
            pstdd[i] = std(φ.prior)
        end

        trspec[i, 1] = φ.transformtype
        (left, right) = φ.transformbounds
        trspec[i, 2] = left
        trspec[i, 3] = right
        if φ == Θ.modelalp_ind
            trspec[i, 4] = 0
        else
            trspec[i, 4] = 1
        end

        i += 1
    end


    para_matlab = readcsv("m990/parameters/para.csv")
    println("### para")
    @test test_matrix_eq(para_matlab, para)

    bounds_matlab = readcsv("m990/parameters/bounds.csv")
    println("### bounds")
    @test test_matrix_eq(bounds_matlab, bounds)

    pshape_matlab = readcsv("m990/parameters/pshape.csv")
    println("### pshape")
    @test test_matrix_eq(pshape_matlab, pshape)

    pmean_matlab = readcsv("m990/parameters/pmean.csv")
    println("### pmean")
    @test test_matrix_eq(pmean_matlab, pmean)
    
    pstdd_matlab = readcsv("m990/parameters/pstdd.csv")
    println("### pstdd")
    @test test_matrix_eq(pstdd_matlab, pstdd)
    
    trspec_matlab = readcsv("m990/parameters/trspec.csv")
    println("### trspec")
    @test test_matrix_eq(trspec_matlab, trspec)
    
    println("parameters.jl tests passed\n")
end



# src/models/m990/modelinds.jl
function test_modelinds()
    # ModelInds object creation
    I = ModelInds(spec_vars())
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
    eq = I.eqconds
    @test length(eq) == 66
    @test eq["eq_Ez"] == 60

    # Additional states
    endo_addl = I.endostates_postgensys
    @test length(endo_addl) == 12
    @test endo_addl["y_t1"] == 67

    # Observables
    obs = I.observables
    @test length(obs) == 18
    @test obs["tfp"] == 12
    
    println("modelinds.jl tests passed\n")
end




# src/models/m990/eqcond.jl
function test_eqcond()
    # eqcond function executes successfully
    model = Model()
    Γ0, Γ1, C, Ψ, Π = eqcond(model)

    # Matrices are of expected dimensions
    @test size(Γ0) == (66, 66)
    @test size(Γ1) == (66, 66)
    @test size(C) == (66, 1)
    @test size(Ψ) == (66, 22)
    @test size(Π) == (66, 13)

    # Check output matrices against Matlab output (ε = 1e-4)
    Γ0_matlab = readcsv("m990/eqcond/Γ0.csv")
    println("### Γ0")
    @test test_matrix_eq(Γ0_matlab, Γ0)

    Γ1_matlab = readcsv("m990/eqcond/Γ1.csv")
    println("### Γ1")
    @test test_matrix_eq(Γ1_matlab, Γ1)

    C_matlab = readcsv("m990/eqcond/C.csv")
    println("### C")
    @test test_matrix_eq(C_matlab, C)

    Ψ_matlab = readcsv("m990/eqcond/PSI.csv")
    println("### Ψ")
    @test test_matrix_eq(Ψ_matlab, Ψ)

    Π_matlab = readcsv("m990/eqcond/PIE.csv")
    println("### Π")
    @test test_matrix_eq(Π_matlab, Π)

    println("eqcond.jl tests passed\n")
end



# src/models/m990/measurement.jl
function test_measurement()
    # TODO: implement

    println("measurement.jl tests passed\n")
end
