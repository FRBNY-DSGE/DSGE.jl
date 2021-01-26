m = PoolModel("ss1")
fp = dirname(@__FILE__)

m <= Setting(:dataroot, "$(fp)/../reference/")
m <= Setting(:data_vintage, "190822")
observables = OrderedDict{Symbol,Observable}()
identity_transform_805 = function(levels)
    return levels[!,:p805]
end
identity_transform_904 = function(levels)
    return levels[!,:p904]
end
observables[:Model904] = Observable(:Model904, [:p904__wrongorigmatlab],
                                    identity_transform_904,
                                    identity_transform_904,
                                    "Model 904 predictive density",
                                    "Model 904 conditional predictive density scores")
observables[:Model805] = Observable(:Model805, [:p805__wrongorigmatlab],
                                    identity_transform_805,
                                    identity_transform_805,
                                    "Model 805 predictive density",
                                    "Model 805 conditional predictive density scores")

@testset "Test various forms of data loading" begin
    m.observable_mappings = observables
    @info "The following warnings are expected"
    df1 = load_data(m; try_disk = false)
    df1 = load_data(m)
    df2 = CSV.read(joinpath(dataroot(m), "raw/wrongorigmatlab_190822.csv"), DataFrame)
    df2[!,:date] = Vector{Dates.Date}(df2[!,:date])
    df2[!,:p904] = Vector{Float64}(df2[!,:p904])
    df2[!,:p805] = Vector{Float64}(df2[!,:p805])

    @test df1[!,:Model904] == df2[!,:p904]
    @test df1[!,:Model805] == df2[!,:p805]
end

nothing
