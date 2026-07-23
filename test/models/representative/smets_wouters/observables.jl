sw = SmetsWouters()

obs = sw.observable_mappings

@testset "observables: expected set of keys" begin
    expected_keys = [:obs_gdp, :obs_hours, :obs_wages, :obs_gdpdeflator,
                     :obs_nominalrate, :obs_consumption, :obs_investment]
    @test collect(keys(obs)) == expected_keys
    @test length(obs) == 7
end

@testset "observables: index mapping" begin
    @test sw.observables[:obs_gdp]          == 1
    @test sw.observables[:obs_hours]        == 2
    @test sw.observables[:obs_wages]        == 3
    @test sw.observables[:obs_gdpdeflator]  == 4
    @test sw.observables[:obs_nominalrate]  == 5
    @test sw.observables[:obs_consumption]  == 6
    @test sw.observables[:obs_investment]   == 7
end

@testset "observables: input series" begin
    @test :GDP__FRED      in obs[:obs_gdp].input_series
    @test :GDPDEF__FRED   in obs[:obs_gdp].input_series

    @test :AWHNONAG__FRED in obs[:obs_hours].input_series
    @test :CE16OV__FRED   in obs[:obs_hours].input_series

    @test :COMPNFB__FRED  in obs[:obs_wages].input_series
    @test :GDPDEF__FRED   in obs[:obs_wages].input_series

    @test :GDPDEF__FRED   in obs[:obs_gdpdeflator].input_series

    @test :DFF__FRED      in obs[:obs_nominalrate].input_series

    @test :PCE__FRED      in obs[:obs_consumption].input_series

    @test :FPI__FRED      in obs[:obs_investment].input_series
end

@testset "observables: names" begin
    @test obs[:obs_gdp].name         == "Real GDP Growth"
    @test obs[:obs_hours].name       == "Hours Per Capita"
    @test obs[:obs_wages].name       == "Percent Change in Wages"
    @test obs[:obs_gdpdeflator].name == "GDP Deflator"
    @test obs[:obs_nominalrate].name == "Nominal FFR"
    @test obs[:obs_consumption].name == "Consumption growth per capita"
    @test obs[:obs_investment].name  == "Real Investment per capita"
end

@testset "observables: rev_transforms" begin
    @test obs[:obs_gdp].rev_transform        == loggrowthtopct_annualized_percapita
    @test obs[:obs_hours].rev_transform       == logleveltopct_annualized_percapita
    @test obs[:obs_wages].rev_transform       == loggrowthtopct_annualized
    @test obs[:obs_gdpdeflator].rev_transform == loggrowthtopct_annualized
    @test obs[:obs_nominalrate].rev_transform == quartertoannual
    @test obs[:obs_consumption].rev_transform == loggrowthtopct_annualized_percapita
    @test obs[:obs_investment].rev_transform  == loggrowthtopct_annualized_percapita
end

nothing
