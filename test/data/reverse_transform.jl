# Load data to use for tests
path = dirname(@__FILE__)
fred = CSV.read("$path/../reference/fred_160812.csv")

# Specify vintage and dates
custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :n_anticipated_shocks     => Setting(:n_anticipated_shocks, 6),
    :population_forecast      => Setting(:population_forecast, false))
m = Model990(testing = true, custom_settings = custom_settings)
m <= Setting(:rate_expectations_source, :ois)

# Test full, semi, and none conditional data
exp_data, exp_cond_data, exp_semicond_data =
    JLD2.jldopen("$path/../reference/load_data_out.jld2", "r") do file
        read(file, "data"), read(file, "cond_data"), read(file, "semi_cond_data")
    end
exp_data_rev_transforms =
    JLD2.jldopen("$path/../reference/transform_data_out.jld2", "r") do file
        read(file, "rev_transform")
    end
pop_forecast = DataFrame(:CNP16OV = cumsum(fred[:CNP16OV][end] .* (ones(8) .* 1.03)))
pop_forecast[:date] = DSGE.get_quarter_ends(DSGE.nextquarter(fred[:date][end]), fred[:date][end] + Year(2))

histobs = load(mbhistobs)
histobs_mat = Matrix{Float64}(histobs[setdiff(names(histobs), :date)])
histpseudo = load(mbhistpseudo)
histpseudo_mat = Matrix{Float64}(histpseudo[setdiff(names(histpseudo), :date)])
reverse_transform(m, histobs_mat, histobs[:date][1], [:obs_gdp, :obs_corepce], :obs)
reverse_transform(m, histpseudo_mat, histpseudo[:date][1], [:OutputGap, :Inflation], :pseudo)
reverse_transform(m, histobs, :obs)
reverse_transform(m, histpseudo, :pseudo)
reverse_transform(m, histobs_mat, histobs[:date][1], [:obs_gdp, :obs_corepce], :obs, fourquarter = true)
reverse_transform(m, histpseudo_mat, histpseudo[:date][1], [:OutputGap, :Inflation], :pseudo,
                  fourquarter = true)
reverse_transform(m, histobs, :obs, fourquarter = true)
reverse_transform(m, histpseudo, :pseudo, fourquarter = true)

transform_list = [loggrowthtopct_4q_percapita, loggrowthtopct_4q, loggrowth_4q_approx,
                  logleveltopct_4q_percapita, logleveltopct_4q, logleveltopct_4q_approx,
                  quartertoannaul, identity]
transform_list_out = load("$path/../reference/transform_data_out,jld2", "transform_list_out")
for (i,transform_fcn) in enumerate(transform_list)
    @test transform_list_out[i] == reverse_transform(histobs_mat[1], transform_fcn, fourquarter = false) # already tested y0s vs. y0
    @test transform_list_out[i] == reverse_transform(histobs_mat[1], transform_fcn, fourquarter = true)
end
@test_throws ErrorException reverse_transform(histobs_mat[1], sum, fourquarter = false)
@test_throws ErrorException reverse_transform(histobs_mat[1], sum, fourquarter = true)
rowth, replace with zeros, test all functions

# Test reverse_transform when passed input and cond type all at once
