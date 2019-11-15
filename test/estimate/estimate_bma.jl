# Set up PoolModel
pm = PoolModel("ss1")
pm <= Setting(:data_vintage, "190822")
filepath = dirname(@__FILE__)
pm <= Setting(:dataroot, "$(filepath)/../reference/")
df = load_data(pm)
data = df_to_matrix(pm, df)

# Run test
bma_ans = zeros(2)
bma_ans[1] = 0.5 * data[1,1] / (.5 * data[1,1] + .5 * data[2,1])
bma_ans[2] = bma_ans[1] * data[1,2] / (bma_ans[1] * data[1,2] + (1 - bma_ans[1]) * data[2,2])
@testset "Check that BMA is correct for PoolModels" begin
    estimate_bma(pm, df[1:2,:]; save_output = false)
    global λ, ~ = estimate_bma(pm, df[1:2,:]; save_output = false, return_output = true)
    @test λ == bma_ans
    global λ, ~ = estimate_bma(pm, data[:,1:2]; save_output = false, return_output = true)
    @test λ == bma_ans
end
