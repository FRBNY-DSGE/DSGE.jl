# To be removed after running this test individually in the REPL successfully
using DSGE
using HDF5, JLD
import Base.Test: @test, @testset

path = dirname(@__FILE__)

m = AnSchorfheide()

n_free_para = count([!m.parameters[i].fixed for i in 1:n_parameters(m)])
free_para_inds = find(x -> x.fixed == false, m.parameters)
n_blocks = 3

srand(42)
test_blocks_free = DSGE.generate_free_blocks(n_free_para, n_blocks)
test_blocks_all  = DSGE.generate_all_blocks(test_blocks_free, free_para_inds)

saved_blocks_free = load("../../reference/util.jld", "blocks_free")
saved_blocks_all  = load("../../reference/util.jld", "blocks_all")

@testset "Mutation block generation" begin
    @test test_blocks_free == saved_blocks_free
    @test test_blocks_all  == saved_blocks_all
end
