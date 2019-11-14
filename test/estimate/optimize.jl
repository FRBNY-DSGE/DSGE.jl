using DSGE
using HDF5
path = dirname(@__FILE__)
writing_output = false

#-----------------------------------------------------------------
# csminwel
#-----------------------------------------------------------------
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

# Load data
file = "$path/../reference/optimize.h5"
x0   = h5read(file, "params")
data = h5read(file, "data")'

# For regenerating test file
params_test = deepcopy(x0)
data_test   = Matrix{Float64}(data)

minimizer  = h5read(file, "minimizer")
minimum    = h5read(file, "minimum")
H_expected = h5read(file, "H")

# See src/estimate/estimate.jl
DSGE.update!(m, x0)
n_iterations = 3

x0 = Float64[p.value for p in m.parameters]
out, H = optimize!(m, data; iterations=n_iterations)

# Re-generate test file
#=if writing_output
    h5open("$path/../reference/optimize.h5", "w") do file
        file["params"] = params_test
        file["data"] = data_test
        file["minimizer"] = out.minimizer
        file["minimum"] = out.minimum
        file["H"] = H
    end
end=#
@testset "Check optimize minimizers are the same [csminwel]" begin
    @test @test_matrix_approx_eq minimizer out.minimizer
    @test minimum ≈ out.minimum atol=1e-6
    @test @test_matrix_approx_eq H_expected H
end

#-----------------------------------------------------------------
# Simulated Annealing
#-----------------------------------------------------------------
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

# Load data
file = "$path/../reference/optimize.h5"#_simulated_annealing.h5"
x0   = h5read(file, "params")
data = h5read(file, "data")'

# For regenerating test file
params_test = deepcopy(x0)
data_test   = Matrix{Float64}(data)

minimizer  = h5read(file, "minimizer")
minimum    = h5read(file, "minimum")
H_expected = h5read(file, "H")

# See src/estimate/estimate.jl
DSGE.update!(m, x0)
n_iterations = 1000

x0 = Float64[p.value for p in m.parameters]
out, H = optimize!(m, data; iterations = n_iterations, method = :simulated_annealing)

# Re-generate test file
if writing_output
    h5open("$path/../reference/optimize_simulated_annealing.h5", "w") do file
        file["params"] = params_test
        file["data"] = data_test
        file["minimizer"] = out.minimizer
        file["minimum"] = out.minimum
        file["H"] = H
    end
end
@testset "Check optimize minimizers are the same [simulated_annealing]" begin
    #@test @test_matrix_approx_eq minimizer out.minimizer
    #@test minimum ≈ out.minimum atol=1e-6
    #@test @test_matrix_approx_eq H_expected H
end

#-----------------------------------------------------------------
# Nelder Mead
#-----------------------------------------------------------------
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

# Load data
file = "$path/../reference/optimize_nelder_mead.jld2"
x0   = load(file, "params")
data = load(file, "data")

# For regenerating test file
params_test = deepcopy(x0)
data_test   = data

minimizer  = load(file, "minimizer")
minimum    = load(file, "minimum")
H_expected = load(file, "H")

# See src/estimate/estimate.jl
DSGE.update!(m, x0)
n_iterations = 3

x0 = Float64[p.value for p in m.parameters]
out, H = optimize!(m, data; iterations = n_iterations, method = :nelder_mead)

# Re-generate test file
if writing_output
    JLD2.jldopen("$path/../reference/optimize_nelder_mead.jld2",
                 true, true, true, IOStream) do file
        file["params"] = params_test
        file["data"] = data_test
        file["minimizer"] = out.minimizer
        file["minimum"] = out.minimum
        file["H"] = H
    end
end
@testset "Check optimize minimizers are the same [nelder_mead]" begin
    @test @test_matrix_approx_eq minimizer out.minimizer
    @test minimum ≈ out.minimum atol=1e-6
    @test @test_matrix_approx_eq H_expected H
end


#-----------------------------------------------------------------
# LBFGS
#-----------------------------------------------------------------
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

# Load data
file = "$path/../reference/optimize.h5"#_lbfgs.h5"
x0   = h5read(file, "params")
data = h5read(file, "data")'

# For regenerating test file
params_test = deepcopy(x0)
data_test   = Matrix{Float64}(data)

minimizer  = h5read(file, "minimizer")
minimum    = h5read(file, "minimum")
H_expected = h5read(file, "H")

# See src/estimate/estimate.jl
DSGE.update!(m, x0)
n_iterations = 3

x0 = Float64[p.value for p in m.parameters]
#out, H = optimize!(m, data; iterations = n_iterations, method = :lbfgs)

# Re-generate test file
if writing_output
    h5open("$path/../reference/optimize_lbfgs.h5", "w") do file
        file["params"] = params_test
        file["data"] = data_test
        file["minimizer"] = out.minimizer
        file["minimum"] = out.minimum
        file["H"] = H
    end
end
@testset "Check optimize minimizers are the same [lbfgs]" begin
    #@test @test_matrix_approx_eq minimizer out.minimizer
    #@test minimum ≈ out.minimum atol=1e-6
    #@test @test_matrix_approx_eq H_expected H
end


nothing
