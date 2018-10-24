using JLD
using Plots

x = load("ct_results.jld")
max_σ = x["max_σ"]
σ_vec1 = x["σ_vec"]
x = load("sim_results.jld")
max_σ_mat = x["max_σ_mat"]
σ_vec2 = x["σ_vec"]

histogram(σ_vec1, max_σ, bins = σ_vec1, title = "ODE Tsitouras", xlabel = "sigma", ylabel = "Freq")
histogram(σ_vec2, max_σ_mat[:, 1], bins = σ_vec2, title = "2 Subintervals", xlabel = "sigma", ylabel = "Freq")
histogram(σ_vec2, max_σ_mat[:, 2], bins = σ_vec2, title = "3 Subintervals", xlabel = "sigma", ylabel = "Freq")
histogram(σ_vec2, max_σ_mat[:, 3], bins = σ_vec2, title = "12 Subintervals", xlabel = "sigma", ylabel = "Freq")
