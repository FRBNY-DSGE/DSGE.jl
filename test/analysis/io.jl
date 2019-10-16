path = dirname(@__FILE__)
using DSGE

m = AnSchorfheide()

@testset "Testing getting MeansBands Input and Output Files" begin
    @test get_meansbands_input_file(m, :mode, :none, :histobs) == rawpath(m, "forecast")*"/histobs_cond=none_para=mode_vint="*data_vintage(m)*".jld2"
    @test get_meansbands_input_file("a", ["b"], :mode, :none, :histobs, fileformat = :jld2) == "a/histobs_b_cond=none_para=mode.jld2"

    @test get_meansbands_output_file(m, :mode, :none, :histobs) == workpath(m, "forecast")*"/mbhistobs_cond=none_para=mode_vint="*data_vintage(m)*".jld2"
    @test get_meansbands_output_file("a", ["b"], :mode, :none, :histobs, fileformat = "jld2") == "a/mbhistobs_b_cond=none_para=mode.jld2"

    @test Matrix(read_mb("$path/../reference/mbhistobs.jld2").means) == Matrix(load("$path/../reference/mbhistobs.jld2", "mb").means)
    @test_throws AssertionError read_bdd_and_unbdd_mb("", "")
    @test Matrix(read_bdd_and_unbdd_mb("$path/../reference/mbbddforecastobs.jld2", "$path/../reference/mbforecastobs.jld2").means) == Matrix(load("$path/../reference/bdd_unbdd.jld2", "bdd_unbdd"))
end

@testset "Testing add_requisite_output_vars_meansbands" begin
    @test DSGE.add_requisite_output_vars_meansbands([:histobs, :shockdecpseudo]) == [:histobs, :shockdecpseudo, :dettrendpseudo, :trendpseudo, :histforecastpseudo]
    @test DSGE.add_requisite_output_vars_meansbands([:histobs, :shockdecobs]) == [:histobs, :shockdecobs, :dettrendobs, :trendobs, :histforecastobs]
end
