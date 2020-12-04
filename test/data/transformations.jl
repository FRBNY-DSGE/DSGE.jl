# using DSGE, ModelConstructors, Dates, CSV, FileIO, Random, JLD2, Test, Nullables
# Load data to use for tests and prep model objet
path = dirname(@__FILE__)
data = JLD2.jldopen("$path/../reference/load_data_out.jld2", "r") do file
    read(file, "data")
end
fred = CSV.read("$path/../reference/fred_160812.csv", DataFrame)
custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :n_anticipated_shocks     => Setting(:n_anticipated_shocks, 6))

m = Model990(custom_settings = custom_settings, testing = true)
m <= Setting(:rate_expectations_source, :ois)
hpfilter_out = JLD2.jldopen("$path/../reference/transformations_out.jld2", "r") do file
    read(file, "hpfilter_out")
end
fred[!,:filtered_population] = hpfilter_out
@testset "Check quarter and percapita helper functions" begin
    @test annualtoquarter(4.) == 1.
    @test quartertoannual(1.) == 4.
    @test quartertoannualpercent(1.) == 400.
    @test @test_matrix_approx_eq (nominal_to_real(:GDP, fred, deflator_mnemonic = :GDPDEF)) (fred[!,:GDP] ./ fred[!,:GDPDEF])
    @test (percapita(:GDP, fred, :CNP16OV)) ≈ (fred[!,:GDP] ./ fred[!,:CNP16OV])
    @test @test_matrix_approx_eq (percapita(m, :GDP, fred))  (fred[!,:GDP] ./ fred[!,:filtered_population])
    m <= Setting(:hpfilter_population, false)
    @test @test_matrix_approx_eq (percapita(m, :GDP, fred))  (fred[!,:GDP] ./ fred[!,:CNP16OV])
    m <= Setting(:population_mnemonic, Nullable())
    @test_throws ErrorException percapita(m, :GDP, fred)
    @test DSGE.lag(fred[!,:GDP], 10)[11:end] ==  fred[!,:GDP][1:end-10]
    hpfilter_test, ~ = hpfilter(fred[!,:CNP16OV], 1600.)
    @test @test_matrix_approx_eq hpfilter_test hpfilter_out
    @test @test_matrix_approx_eq difflog([1., 3., 1.]) vec([missing, log(3), -log(3)])
    @test @test_matrix_approx_eq oneqtrpctchange([1., 3., 1.])  100 * vec([missing, log(3), -log(3)])
end


@testset "Check loggrowth and loglevel transforms" begin
    @test loggrowthtopct([1., 3.]) == 100 .* (exp.([.01, .03]) .- 1.)
    @test loggrowthtopct_percapita([1., 3.], [1., 3.]) == 100. .* ((exp.([.01, .03]) .*
                                                                   exp.([1., 3.]).^4) .- 1.)
    @test loggrowthtopct_percapita([1. 3.; 1. 3.], [1., 3.]) == (100. .* ((exp.([.01, .03]) .*
                                                                           exp.([1., 3.]).^4) .- 1.))' .* ones(2)
    @test_throws AssertionError loggrowthtopct_percapita([1., 3.], [1., 3., 4.])
    @test loggrowthtopct_annualized([1., 3.]) == 100. .* (exp.([.01, .03]).^4 .- 1.)
    @test loggrowthtopct_annualized_percapita([2., 3.], [1.1, 1.2]) == 100. * (exp.([2., 3.] ./ 100. .+  [1.1, 1.2]).^4 .- 1.)
    @test loggrowthtopct_annualized_percapita([2. 3.; 2. 3.], [1.1, 1.2]) == ones(2) * (100. * (exp.([2., 3.] ./ 100. .+
                                                                                                  [1.1, 1.2]).^4 .- 1.))'
    @test_throws AssertionError loggrowthtopct_annualized_percapita([1., 3.], [1., 3., 4.])

    @test isnan(logleveltopct_annualized_approx([.1, .2])[1])
    @test logleveltopct_annualized_approx([.1, .2], 0.) == [.4, .4]
    @test sum(isnan.(logleveltopct_annualized_approx([.1 .2; 1. 2.]))) == 2
    @test logleveltopct_annualized_approx([.1 .2; 1. 2.], 0.) == [.4 .4; 4 4]
end

@testset "Check retrieval functions for transforms" begin
    @test DSGE.get_nopop_transform(loggrowthtopct_annualized_percapita) == loggrowthtopct_annualized
    @test DSGE.get_nopop_transform(logleveltopct_annualized_percapita) == logleveltopct_annualized
    @test DSGE.get_nopop_transform(sum) == sum
    @test DSGE.get_irf_transform(loggrowthtopct_annualized_percapita) == loggrowthtopct_annualized
    @test DSGE.get_irf_transform(logleveltopct_annualized_percapita) == logleveltopct_annualized
    @test DSGE.get_irf_transform(sum) == sum
    @test DSGE.get_transform4q(loggrowthtopct_annualized_percapita) == DSGE.loggrowthtopct_4q_percapita
    @test DSGE.get_transform4q(loggrowthtopct_annualized) == DSGE.loggrowthtopct_4q
    @test DSGE.get_transform4q(logleveltopct_annualized_percapita) == DSGE.logleveltopct_4q_percapita
    @test DSGE.get_transform4q(logleveltopct_annualized) == DSGE.logleveltopct_4q
    @test DSGE.get_transform4q(identity) == identity
    @test_throws ErrorException DSGE.get_transform4q(sum)

    @test DSGE.get_scenario_transform(loggrowthtopct_annualized_percapita) == quartertoannual
    @test DSGE.get_scenario_transform(loggrowthtopct_annualized) == quartertoannual
    @test DSGE.get_scenario_transform(quartertoannual) == quartertoannual
    @test DSGE.get_scenario_transform(logleveltopct_annualized_percapita) == logleveltopct_annualized_approx
    @test DSGE.get_scenario_transform(logleveltopct_annualized) == logleveltopct_annualized_approx
    @test DSGE.get_scenario_transform(identity) == identity
    @test_throws ErrorException DSGE.get_scenario_transform(sum)
end

@testset "Check 4q transforms" begin
    @test DSGE.prepend_data([1., 2.], [1., 2.]) == [1., 2., 1., 2.]
    @test DSGE.prepend_data([1. 2.; 1. 2.], [1.; 2. ]) == [1. 2. 1. 2. ; 1. 2. 1. 2.]
    @test sum(isnan.(DSGE.loggrowthtopct_4q([1., 2.]))) == 2
    @test DSGE.loggrowthtopct_4q([1., 2.], [1., 1., 1.]) == 100. .* (exp.([4., 5.] ./ 100) .- 1.)
    @test sum(isnan.(DSGE.loggrowthtopct_4q([1. 2.; 1. 2.]))) == 4
    @test DSGE.loggrowthtopct_4q([1. 2.; 1. 2.], [1., 1., 1.]) == 100. .* (exp.([4. 5.; 4. 5.] ./ 100) .- 1.)
    @test sum(isnan.(DSGE.loggrowthtopct_4q_percapita([1.], [1., 1., 1., 1.]))) == 1
    @test DSGE.loggrowthtopct_4q_percapita([1.], [1., 1., 1., 1.], [1., 1., 1.]) == 100. .* (exp.([4.] ./ 100. .+ [4.]) .- 1.)
    @test sum(isnan.(DSGE.loggrowthtopct_4q_percapita([1. 1.; 1. 1.], [1., 1., 1., 1., 1.]))) == 4
    @test DSGE.loggrowthtopct_4q_percapita([1. 1.; 1. 1.], [1., 1., 1., 1., 1.], [1., 1., 1.]) == ones(2) * 100. .* (exp.([4.; 4.] ./
                                                                                                                100. .+ [4.; 4.]) .- 1.)'
    @test sum(isnan.(DSGE.logleveltopct_4q([1., 2., 3., 4.]))) == 4
    @test DSGE.logleveltopct_4q([1., 2., 3., 4. ], [1., 1., 1., 1.]) == 100. .* (exp.([0., 1., 2., 3.] ./ 100) .- 1.)
    @test sum(isnan.(DSGE.logleveltopct_4q([1. 2. 3. 4. ; 1. 2. 3. 4. ]))) == 8
    @test DSGE.logleveltopct_4q([1. 2. 3. 4. ; 1. 2. 3. 4.], [1., 1., 1., 1.]) == 100. .* (exp.([0. 1. 2. 3.; 0. 1. 2. 3.] ./ 100) .- 1.)
    @test sum(isnan.(DSGE.logleveltopct_4q_percapita([1., 2., 3., 4.], ones(7)))) == 4
    @test DSGE.logleveltopct_4q_percapita([1., 2., 3., 4.], ones(7), zeros(4)) == 100. .* (exp.([1., 2.,
                                                                                                 3., 4.] ./ 100. .+
                                                                                                [4., 4., 4., 4.]) .- 1.)
    @test sum(isnan.(DSGE.logleveltopct_4q_percapita([1. 2. 3. 4.; 1. 2. 3. 4.], ones(7)))) == 8
    @test DSGE.logleveltopct_4q_percapita([1. 2. 3. 4.; 1. 2. 3. 4.], ones(7), zeros(4)) == 100. .*
              (exp.([1. 2. 3. 4.; 1. 2. 3. 4.] ./ 100. + 4 * ones(2,4)) .- 1.)
    @test_throws AssertionError DSGE.logleveltopct_4q_percapita([1., 2., 3., 4.], [1., 2., 3.])

    @test sum(isnan.(loggrowthtopct_4q_approx([.1, .2]))) == 2
    @test @test_matrix_approx_eq loggrowthtopct_4q_approx([.1, .2], zeros(3)) [.1, .3]
    @test @test_matrix_approx_eq loggrowthtopct_4q_approx([.1 .2; .1 .2], zeros(3)) [.1 .3; .1 .3]
    @test_throws AssertionError loggrowthtopct_4q_approx([.1, .2], [.1, .2])

    @test sum(isnan.(DSGE.logleveltopct_4q_approx([.1, .2, .3, .4]))) == 4
    @test @test_matrix_approx_eq DSGE.logleveltopct_4q_approx([.1, .2, .3, .4], [.1, .1, .1, .1]) [0., .1, .2, .3]
    @test_throws AssertionError DSGE.logleveltopct_4q_approx([.1, .2], [.1, .2])
end

nothing
