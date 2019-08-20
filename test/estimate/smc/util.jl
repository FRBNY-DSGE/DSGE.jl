path = dirname(@__FILE__)

m = AnSchorfheide()

n_free_para = count([!m.parameters[i].fixed for i in 1:n_parameters(m)])
free_para_inds = findall(x -> x.fixed == false, m.parameters)
n_blocks = 3

@everywhere Random.seed!(42)
test_blocks_free = SMC.generate_free_blocks(n_free_para, n_blocks)
test_blocks_all  = SMC.generate_all_blocks(test_blocks_free, free_para_inds)

saved_blocks_free = load("$path/../../reference/util.jld2", "blocks_free")
saved_blocks_all  = load("$path/../../reference/util.jld2", "blocks_all")

@testset "Mutation block generation" begin
    @test test_blocks_free == saved_blocks_free
    @test test_blocks_all  == saved_blocks_all
end
