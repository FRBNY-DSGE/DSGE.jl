# Load data to use for tests and prep model objet
path = dirname(@__FILE__)
data = JLD2.jldopen("$path/../reference/load_data_out.jld2", "r") do file
    read(file, "data")
end
fred = CSV.read("$path/../reference/fred_160812.csv")
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

@testset "Check quarter and percapita helper functions" begin
    @test annualtoquarter(4.) == 1.
    @test quartertoannual(1.) == 4.
    @test quatertoannualpercent(1.) == 400.
    @test (nominal_to_real(:GDP, fred, deflator_mnemonic = :GDPDEF)) ≈ (fred[:GDP] ./ fred[:GDPDEF])
    @test (percapita(:GDP, fred, population_mnemonic = :CNP16OV)) ≈ s(fred[:GDP] ./ fred[:CNP16OV])
    @test (percapita(m, :GDP, fred)) ≈ (fred[:GDP] ./ fred[:CNP16OV])
    m <= Setting(:hpfilter_population, false)
    @test (percapita(m, :GDP, fred)) ≈ (fred[:GDP] ./ fred[:CNP16OV])
    m <= Setting(:population_mnemonic, Nullable())
    @test_throws ErrorException percapita(m, :GDP, fred)
    @test lag(fred[:GDP], 10) == vcat(NaN*ones(10), fred[:GDP][1:end-10])
    hpfilter_out = JLD2.jldopen("$path/../reference/transformations_out.jld2", "r") do file
        read(file, "hpfilter_out")
    end
    @test hpfilter(fred[:GDP], 1600.) ≈ hpfilter_out
    @test difflog([1., 3., 1.]) == vec([missing, log(3), -log(3)])
    @test difflog(hcat([1., 3., 1.], [1., 1., 1.])) == vec([missing log(3) -log(3) 0. 0. 0.])
    @test oneqtrpctchange([1., 3., 1.]) = 100 * vec([missing, log(3), -log(3)])
end

@testset "Check loggrowth and loglevel transforms" begin
    @test loggrowthtopct([1., 3.]) = 100 .* (exp.([.01, .03]) .- 1.)
    @test loggrowthtopct_percapita([1., 3.], [1., 3.]) == 100. .* ((exp.([.01, .03]) .*
                                                                   exp.([1., 3.]).^4) - 1.)
    @test loggrowthtopct_percapita([1. 3.; 1. 3.], [1., 3.]) == (100. .* ((exp.([.01, .03]) .*
                                                                           exp.([1., 3.]).^4) - 1.))' .* ones(2)
    @test_throws AssertionError loggrowthtopct_percapita([1., 3.], [1., 3., 4.])
    @test loggrowthtopct_annualized([1., 3.]) == 100. .* (exp.([.01, .03]).^4 .- 1.)
    @test loggrowthtopct_annualized_percapita([2., 3.], [1.1, 1.2], 1.) == 100. * (exp.([2., 3.] ./ 100. -
                                                                                      [1., 2.] ./ 100. .+
                                                                                      [1.1, 1.2]).^4 .- 1.)
    @test loggrowthtopct_annualized_percapita([2., 3.; 2. 3.], [1.1, 1.2], 1.) == ones(2) .* (100. * (exp.([2., 3.] ./ 100. -
                                                                                                  [1., 2.] ./ 100. .+
                                                                                                  [1.1, 1.2]).^4 .- 1.))
    @test loggrowthtopct_annualized_percapita([2., 3.], [1.1, 1.2]) == 100. * (exp.([2., 3.] ./ 100. -
                                                                                  [NaN, 2.] ./ 100. .+
                                                                                  [1.1, 1.2]).^4 .- 1.)
    @test loggrowthtopct_annualized_percapita([2., 3.], [1.1, 1.2]) == ones(2) .* (100. * (exp.([2., 3.] ./ 100. -
                                                                                              [NaN, 2.] ./ 100. .+
                                                                                              [1.1, 1.2]).^4 .- 1.))
    @test_throws AssertionError loggrowthtopct_annualized_percapita([1., 3.], [1., 3., 4.])

    @test logleveltopct_annualized_approx([.1, .2]) == [NaN, .4]
    @test logleveltopct_annualized_approx([.1, .2], 0.) == [.4, .4]
    @test logleveltopct_annualized_approx([.1 .2; 1. 2.]) == [NaN .4; NaN .4]
    @test logleveltopct_annualized_approx([.1 .2; 1. 2.], 0.) == [.4 .4; .4 .4]
end

@testset "Check retrieval functions for transforms" begin
    @test get_nopop_transform(loggrowthtopct_annualized_percapita) == loggrowthtopct_annualized
    @test get_nopop_transform(logleveltopct_annualized_percapita) == logleveltopct_annualized
    @test get_nopop_transform(sum) == sum
    @test get_irf_transform(loggrowthtopct_annualized_percapita) == loggrowthtopct_annualized
    @test get_irf_transform(logleveltopct_annualized_percapita) == logleveltopct_annualized
    @test get_irf_transform(sum) == sum
    @test get_transform4q(loggrowthtopct_annualized_percapita) == loggrowthtopct_4q
    @test get_transform4q(loggrowthtopct_annualized) == loggrowthtopct_4q
    @test get_transform4q(logleveltopct_annualized_percapita) == logleveltopct_4q
    @test get_transform4q(logleveltopct_annualized) == logleveltopct_4q
    @test get_transform4q(identity) == identity
    @test_throws ErrorException get_transform4q(sum)

    @test get_scenario_transform(loggrowthtopct_annualized_percapita) == quartertoannual
    @test get_scenario_transform(loggrowthtopct_annualized) == quartertoannual
    @test get_scenario_transform(quartertoannual) == quartertoannual
    @test get_secnario_transform(logleveltopct_annualized_percapita) == logleveltopct_annualized_approx
    @test get_secnario_transform(logleveltopct_annualized) == logleveltopct_annualized_approx
    @test get_scenario_transform(identity) == identity
    @test_throws ErrorException get_scenario_transform(sum)
end

@testset "Check 4q transforms" begin
    @test prepend_data([1., 2.], [1., 2.]) == [1., 2., 1., 2.]
    @test prepend_data([1. 2.; 1. 2.], [1.; 2. ]) == [1. 2. 1. ; 1. 2. 2.]
    @test loggrowthtopct_4q([1., 2.]) == [NaN, NaN]
    @test loggrowthtopct_4q([1., 2.], [1., 1., 1.]) == 100. .* (exp.([4., 5.] ./ 100) .- 1.)
    @test loggrowthtopct_4q([1. 2.; 1. 2.]) == [NaN NaN; NaN NaN]
    @test loggrowthtopct_4q([1. 2.; 1. 2.], [1., 1., 1.]) == 100. .* (exp.([4. 5.; 4. 5.] ./ 100) .- 1.)
    @test loggrowthtopct_4q_percapita([1., 2.], [1., 1.]) == [NaN, NaN]
    @test loggrowthtopct_4q_percapita([1., 2.], [1., 1.], [1. 1. 1.]) == 100. .* (exp.([4., 5.] ./ 100. .+ [1., 1.]) .- 1.)
    @test loggrowthtopct_4q_percapita([1. 2.; 1. 2.], [1., 1.]) == [NaN NaN; NaN NaN]
    @test loggrowthtopct_4q_percapita([1. 2.; 1. 2.], [1., 1.], [1. 1. 1.]) == 100. .* (exp.([4. 5.; 4. 5.] ./ 100. .+ [1. 1.; 1. 1.]) .- 1.)
    @test logleveltopct_4q([1., 2., 3., 4.]) == [NaN, NaN]
    @test logleveltopct_4q([1., 2., 3., 4. ], [1., 1., 1., 1.]) == 100. .* (exp.([0., 1., 2., 3.] ./ 100) .- 1.)
    @test logleveltopct_4q([1. 2. 3. 4. ; 1. 2. 3. 4.]) == [NaN NaN; NaN NaN]
    @test logleveltopct_4q([1. 2. 3. 4. ; 1. 2. 3. 4.], [1., 1., 1.]) == 100. .* (exp.([0. 1. 2. 3.; 0. 1. 2. 3.] ./ 100) .- 1.)
    @test logleveltopct_4q_percapita([1., 2., 3., 4.], [1., 1.]) == [NaN, NaN]
    @test logleveltopct_4q_percapita([1., 2., 3., 4.], [1., 1., 3., 4.], zeros(8)) == 100. .* (exp.([0., 1., 2., 3.] ./ 100.) .- 1.)
    @test logleveltopct_4q_percapita([1. 2. 3. 4.; 1. 2. 3. 4.], [1., 1.]) == [NaN NaN; NaN]
    @test logleveltopct_4q_percapita([1. 2. 3. 4.; 1. 2. 3. 4.], [1., 1.], zeros(8)) == 100. .* (exp.([0. 1. 2. 3.; 0. 1. 2. 3.] ./ 100.) .- 1.)
    @test_throws AssertionError logleveltopct_4q_percapita([1., 2., 3., 4.], [1., 2., 3.])

    @test loggrowthtopct_4q_approx([.1, .2]) == [NaN NaN]
    @test loggrowthtopct_4q_approx([.1, .2], zeros(3)) == [.1, .3]
    @test loggrowthtopct_4q_approx([.1 .2; .1 .2], zeros(3)) == [.1 .3; .1 .3]
    @test_throws AssertionError loggrowthtopct_4q_approx([.1, .2], [.1, .2])

    @test logleveltopct_4q_approx([.1, .2, .3, .4]) == [NaN NaN]
    @test logleveltopct_4q_approx([.1, .2, .3, .4], [.1, .1, .1, .1]) == [0., .1, .2, .3]
    @test_throws AssertionError logleveltopct_4q_approx([.1, .2], [.1, .2])
end
