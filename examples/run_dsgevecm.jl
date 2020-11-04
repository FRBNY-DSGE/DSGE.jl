using DSGE, Plots, FileIO
using Plots.PlotMeasures

# This script loads in matrices from the DSGE in Del Negro, Schorfheide, Smets, and Wouters (2007)
# and computes the DSGE impulse responses to exogenous structural shocks and the
# corresponding DSGE-VECM impulse responses with 4 lags using either λ = ∞ or
# the rotation matrix from the DSGE to identify the
# mapping from the innovations in the VECM to the structural shocks.
# The script creates a dictionary called `plots_dict`, whose keys are the names of
# different shocks and values are the plots of the impulse responses of different
# observables to the shock.
#
# Generally, the DSGE-VECM IRFs with λ = 0.33 will look different from the DSGE IRFs, but
# they will yield an overall similar response. See the response of output growth to
# the z shock. In contrast, the λ = ∞ DSGE-VECM IRFs are very close to the DSGE IRFs.
# See the IRF of output growth to the μ or λ_f shocks for examples of DSGE-VECM IRFs
# with λ = ∞ that are nearly identical to the underlying DSGE IRFs. For IRFs that are
# not exactly the same, see the IRF of output growth to the z shock. On impact, the
# IRF is different, but after 5-10 quarters, the impulse responses are almost identical.
fp = dirname(@__FILE__)

matdata     = load(joinpath(fp, "../test/reference/dsgevecm_lambda_irfs.jld2"))
horizon     = 20
nshocks     = size(matdata["RRR"], 2)
obs_dict    = Dict{Symbol, Int}(:output_growth => 1, :consumption_growth => 2,
                                :investment_growth => 3, :hours => 4,
                                :wage_growth => 5, :inflation => 6, :nominalrate => 7)
# The shocks are technology, labor preference, investment price, preference,
# government spending, price markup, and monetary policy
shocks_dict = Dict{Symbol, Int}(:z => 1, :φ => 2, :μ => 3, :b => 4, :g => 5,
                                :λ_f => 6, :r => 7)

# DSGEVECM IRFs in deviations with λ = 0.33
ŷ_dsgevecm_λ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                    vec(matdata["DD"]), matdata["MM"],
                                    matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                    Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                    matdata["cointvec"], horizon, deviations = true)

# DSGEVECM IRFs with λ = ∞, p = 12
β_λ_∞, Σ_λ_∞ = DSGE.vecm_approx_state_space(matdata["TTT"], matdata["RRR"], matdata["QQ"],
                                       vec(matdata["DD"]), matdata["ZZ"],
                                       zeros(Int(matdata["nvar"]) + Int(matdata["coint"]),
                                             Int(matdata["nvar"]) + Int(matdata["coint"])),
                                       matdata["MM"], Int(matdata["nvar"]), 4, Int(matdata["coint"]),
                                       use_intercept = true)
Σ_λ_∞ = (Σ_λ_∞ + Σ_λ_∞') ./ 2
ŷ_dsgevecm_∞ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                      vec(matdata["DD"]), matdata["MM"],
                                      matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                      Int(matdata["coint"]), β_λ_∞, Σ_λ_∞,
                                      matdata["cointvec"], horizon, deviations = true)

# DSGE IRFs in deviations
dsge_states = zeros(size(matdata["TTT"], 1), horizon, nshocks)
ŷ_dsge = zeros(Int(matdata["nvar"]), horizon, nshocks)
for i in 1:nshocks
    shock = zeros(nshocks)
    shock[i] = -sqrt(matdata["QQ"][i, i])
    dsge_states[:, 1, i] = matdata["RRR"] * shock
    for t in 2:horizon
        dsge_states[:, t, i] = matdata["TTT"] * dsge_states[:, t - 1, i]
    end
    ŷ_dsge[:, :, i] = matdata["ZZ"][1:Int(matdata["nvar"]), :] * dsge_states[:, :, i]
end


# Plot!
plots_dict = Dict{Symbol, Dict{Symbol, Plots.Plot}}()
for (shock, j) in  shocks_dict
    plots_dict[shock] = Dict{Symbol, Plots.Plot}()
    for (obs, i) in obs_dict
        plots_dict[shock][obs] = plot(1:horizon, ŷ_dsgevecm_∞[i, :, j],
                                      label = "DSGE-VECM \\lambda = Inf", color = :black,
                                      linewidth = 3)
        plot!(1:horizon, ŷ_dsge[i, :, j], label = "DSGE",
              color = :red, linewidth = 3, linestyle = :dash)
        plot!(1:horizon, ŷ_dsgevecm_λ[i, :, j], label = "DSGE-VECM \\lambda = 0.33", color = :black,
              linestyle = :dash,
              legend = :bottomright, left_margin = 20px)
    end
end
